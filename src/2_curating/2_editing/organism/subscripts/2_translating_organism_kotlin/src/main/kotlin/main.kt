import com.univocity.parsers.tsv.TsvWriter
import com.univocity.parsers.tsv.TsvWriterSettings

import java.io.File
import java.text.DecimalFormat
import java.util.concurrent.atomic.AtomicInteger
import java.util.stream.Collectors
import kotlin.system.exitProcess

// By how much the gnfinder files are going to be cut into

const val cut = 10_000

// Show usage and exit
fun usageExit() {
    println("Usage: java -jar <scriptname> Data_path full/min/test")
    println("If the last argument is test, it will run on tables_test")
    println("If the last argument is min, it will run on tables_min")
    println("If the last argument is full, it will run on tables")
    exitProcess(1)
}

lateinit var exclusionDic: List<Regex>
lateinit var combinedDic: List<Pair<Regex, String>>

inline fun processRecord(record: Map<String, String>): Map<String, String> {
    val newRecord = record.toMutableMap()
    newRecord["organismInterim"] = ""
    var newTerm = newRecord["organismOriginal"] ?: return newRecord

    exclusionDic.map {
        newTerm = newTerm.replace(it, "")
    }

    if (newTerm == "") return newRecord

    val orgCleaned = newRecord["organismCleaned"]
    if (orgCleaned != null) {
        if (orgCleaned == newTerm) newTerm = ""
        newTerm = newTerm.replace(
            Regex("\\b$orgCleaned\\b", RegexOption.IGNORE_CASE),
            ""
        )
    }

    // We try to clean up a bit

    newTerm = newTerm.replace(Regex("[. _]"), " ").replace(Regex(" +"), " ").trim()
    if (newTerm == "") return newRecord
    val candidateList = mutableSetOf<String>()
    // So this is a bit complicated
    // We match each string only once and remove it from the initial string so it doesn't get matched again
    //
    // So for "foo bar bim" if we have in the dictionnary "foo bar" => ABC  and "bim" => XYZ and "foo" => 123
    // it would process it like that: initialString = "foo bar"
    // initialString = "bar" candidatesList=["ABC"]
    // initialString = "" candidatesList=["ABC", "XYZ"]

    // AR comment: Actually, does not work as expected. For some reason "foo bar" is not properly removed from
    // initialString, so foo gets matched and we end up with candidatesList=["ABC", "XYZ", "123"]. I think 
    // it might be because of a case-sensitive issue

    // We need to make a loop we can exit from
    run loop@{
        combinedDic.forEach { pair ->
            // If we have at least one match
            pair.first.find(newTerm)?.let {
                // We remove it
                newTerm.replace(pair.first, "")
                // And add the replacement to the candidateList
                if (pair.second != "")
                    candidateList.add(pair.second)
                newTerm = newTerm.trim() // We remove eventual unnecessary space

                if (newTerm == "") return@loop
            }
        }
    }

    // We likely don't need the second round as we are doing the cleaning ourselves in the dictionnary creation.
    /*// Second round of TCM

    tcmNamesDic.forEach {
    newTerm = newTerm.replace(it.first, it.second)
    }*/
    // We combine both the new term and the eventual candidateList
    newRecord["organismInterim"] = candidateList.joinToString(" , ")
    return newRecord
}

fun main(args: Array<String>) {
    val logger = MyDirtyLogger() //LoggerFactory.getLogger("main")
    if (args.size < 2) usageExit()
    if (args[1] !in setOf("test","min", "full")) usageExit()

    val pathData = args[0]

    val pathDataInterim = "$pathData/interim"

    val pathDataInterimDictionaries = "$pathDataInterim/dictionaries"
    val pathDataInterimDictionariesTcmNames = "$pathDataInterimDictionaries/tcm/names.tsv.gz"
    val pathDataInterimDictionariesCommonNames = "$pathDataInterimDictionaries/common/names.tsv.gz"
    val pathDataInterimDictionariesCommonBlackDic = "$pathDataInterimDictionaries/common/black.tsv"

    val pathDataInterimTables = when (args[1]) {
        "full" -> "$pathData/interim/tables"
        "min" -> "$pathData/interim/tables_min"
        "test" -> "$pathData/interim/tables_test"
        else -> throw Exception("This shouldn't have happened, we only know these types of tables")
    }

    val pathDataInterimTablesCleaned = "$pathDataInterimTables/2_cleaned"
    val pathDataInterimTablesCleanedOrganism = "$pathDataInterimTablesCleaned/organism"
    val pathDataInterimTablesCleanedOrganismOriginalUniqueTable = "$pathDataInterimTablesCleanedOrganism/originalUnique.tsv.gz"
    val pathDataInterimTablesCleanedOrganismTranslatedInterim =
        "$pathDataInterimTablesCleanedOrganism/interim.tsv.gz"

    val pathDataInterimTablesTranslatedOrganism = "$pathDataInterimTables/1_translated/organism"

    logger.info("Making sure directories exist")
    File(pathDataInterimTablesTranslatedOrganism).mkdirs()
    File(pathDataInterimTablesCleanedOrganism).mkdirs()

    // Organism List

    logger.info("Loading and processing the organism list")
    val dataCleanedOriginalOrganism = parseTSVZFile(pathDataInterimTablesCleanedOrganismOriginalUniqueTable)?.map {
        it.toFieldMap()
    } ?: throw Exception("Sorry can't read organism list.")

    // TCM

    logger.info("Loading and processing the TCM names")

    val tcmNamesDic = parseTSVZFile(pathDataInterimDictionariesTcmNames)?.flatMap {
        val vernacularName = it.getValue("vernacularName", "")
        val canonicalName = it.getValue("canonicalName", "")
        val newCanonicalName = it.getValue("newCanonicalName", "")
        val out = mutableListOf(
            Pair(
                Regex("\\b${vernacularName}\\b", RegexOption.IGNORE_CASE),
                if (newCanonicalName != "NA") {
                    newCanonicalName
                } else {
                    canonicalName
                }
            )  // If we have a new canonical name, we replace it already
        )
        // If we find out a canonicalName and we have a newCanonicalName, we will replace it by that
        if (newCanonicalName != "NA")
            out.add(Pair(Regex("\\b${canonicalName}\\b", RegexOption.IGNORE_CASE), newCanonicalName))
        out
    }?.distinct() ?: throw Exception("Sorry can't read TCM names.")

    // Common Names

    logger.info("Loading and processing common names")
    val commonNamesDic = parseTSVZFile(pathDataInterimDictionariesCommonNames)?.map {
        val vernacularName = it.getValue("vernacularName", "")
        val canonicalName = it.getValue("canonicalName", "")
        Pair(Regex("\\b${vernacularName}\\b", RegexOption.IGNORE_CASE), canonicalName)
    }?.distinct() ?: throw Exception("Sorry can't read common names.")

    // Exclusion list

    logger.info("Loading and processing the exclusion list")
    exclusionDic = parseTSVFile(pathDataInterimDictionariesCommonBlackDic)?.map {
        it.getValue<String>("blackName", null)
    }?.sortedByDescending { it.length }?.map { Regex("\\b$it\\b", RegexOption.IGNORE_CASE) }?.distinct()
        ?: throw Exception("Sorry can't read exclusion list")

    // Processing
    logger.info("Creating combined Dic")
    combinedDic = (commonNamesDic + tcmNamesDic)
    val entriesNumber = dataCleanedOriginalOrganism.size
    logger.info("Processing $entriesNumber entries")
    val progress = AtomicInteger(0)
    val startTime = System.nanoTime()

    val processedRecords =
        dataCleanedOriginalOrganism.parallelStream().map { record ->
            val localProgress = progress.incrementAndGet()
            if (localProgress % 1000 == 0) {
                val ratio = localProgress.toFloat() / entriesNumber
                val timeSpentSeconds = (System.nanoTime() - startTime) / 1_000_000_000
                val rest = entriesNumber - localProgress
                val eta = rest * timeSpentSeconds / localProgress
                logger.info("Processed $localProgress/$entriesNumber ${(100 * ratio)}% ETA: $eta s")
            }
            processRecord(record)
        }.unordered().collect(Collectors.toList())


    logger.info("Done processing")

    // Writing files

    logger.info("Writing Interim file $pathDataInterimTablesCleanedOrganismTranslatedInterim")
    val outputWriter = TsvWriter(GZIPWrite(pathDataInterimTablesCleanedOrganismTranslatedInterim), TsvWriterSettings())
    val headers = processedRecords.firstOrNull()?.keys ?: setOf("organismOriginal	value_min","dbQuality","value_max",
        "rank","ids","taxonomy","organismCleaned","organismDbTaxo","taxonId","organismInterim")
    outputWriter.writeHeaders(headers)
    processedRecords.map { row ->
        if (row != null) outputWriter.writeRow(headers.map { row.getOrDefault(it, "") })
    }
    outputWriter.close()

    logger.info("Writing GNFinder files in $pathDataInterimTablesTranslatedOrganism")
    processedRecords.map { it?.get("organismInterim") ?: "" }.filter { it != "" }
        .chunked(cut).mapIndexed { idx, chunk ->
            val fileName =
                "$pathDataInterimTablesTranslatedOrganism/${DecimalFormat("000000").format(idx * cut + cut)}.tsv"
            val writer = TsvWriter(File(fileName), TsvWriterSettings())
            writer.writeHeaders("organismInterim")
            chunk.map { writer.writeRow(it) }
            writer.close()
        }
    logger.info("Completely done. Go Kotlin \\o/")
}