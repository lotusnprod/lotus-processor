import com.univocity.parsers.tsv.TsvWriter
import com.univocity.parsers.tsv.TsvWriterSettings
import java.io.File
import java.util.concurrent.atomic.AtomicInteger
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


fun main(args: Array<String>) {
    val logger = MyDirtyLogger() //LoggerFactory.getLogger("main")
    if (args.size < 2) usageExit()
    if (args[1] !in setOf("test", "min", "full")) usageExit()

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
    val pathDataInterimTablesCleanedOrganismOriginalUniqueTable =
        "$pathDataInterimTablesCleanedOrganism/originalUnique.tsv.gz"
    val pathDataInterimTablesCleanedOrganismTranslatedInterim =
        "$pathDataInterimTablesCleanedOrganism/interim.tsv.gz"

    val pathDataInterimTablesTranslatedOrganism = "$pathDataInterimTables/1_translated/organism"

    logger.info("Making sure directories exist")
    File(pathDataInterimTablesTranslatedOrganism).mkdirs()
    File(pathDataInterimTablesCleanedOrganism).mkdirs()

    // Organism List
    //
    // We load the data from the organism list and transform it as a field Map

    logger.info("Loading and processing the organism list")
    val dataCleanedOriginalOrganism = parseTSVZFile(pathDataInterimTablesCleanedOrganismOriginalUniqueTable)?.map {
        it.toFieldMap()
    } ?: throw Exception("Sorry can't read organism list.")

    // Now we generate the dictionaries for replacement

    val dics = Dictionaries()

    // TCM
    //
    // Here we are going to generate a list of regular expressions that are going to allow us to match TCM names and
    // replace them.

    logger.info("Loading and processing the TCM names")

    dics.combinedDic.addAll(
        parseTSVZFile(pathDataInterimDictionariesTcmNames)?.flatMap {
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
    )

    // Common Names

    logger.info("Loading and processing common names")
    dics.combinedDic.addAll(
        parseTSVZFile(pathDataInterimDictionariesCommonNames)?.map {
            val vernacularName = it.getValue("vernacularName", "")
            val canonicalName = it.getValue("canonicalName", "")
            Pair(Regex("\\b${vernacularName}\\b", RegexOption.IGNORE_CASE), canonicalName)
        }?.distinct() ?: throw Exception("Sorry can't read common names.")
    )

    // Exclusion list

    logger.info("Loading and processing the exclusion list")
    dics.exclusionDic.addAll(parseTSVFile(pathDataInterimDictionariesCommonBlackDic)?.map {
        it.getValue<String>("blackName", null)
    }?.sortedByDescending { it.length }?.map { Regex("\\b$it\\b", RegexOption.IGNORE_CASE) }?.distinct()
        ?: throw Exception("Sorry can't read exclusion list")
    )

    val entriesNumber = dataCleanedOriginalOrganism.size
    logger.info("Processing $entriesNumber entries")
    val progress = AtomicInteger(0)
    val startTime = System.nanoTime()

    dataCleanedOriginalOrganism.parallelStream().forEach { record ->
        val localProgress = progress.incrementAndGet()
        if (localProgress % 1000 == 0) {
            val ratio = localProgress.toFloat() / entriesNumber
            val timeSpentSeconds = (System.nanoTime() - startTime) / 1_000_000_000
            val rest = entriesNumber - localProgress
            val eta = rest * timeSpentSeconds / localProgress
            logger.info("Processed $localProgress/$entriesNumber ${(100 * ratio)}% ETA: $eta s")
        }

        record["organismInterim"] = dics.processRecord(record)
    }

    logger.info("Done processing the entries")

    // Writing files

    logger.info("Writing Interim file $pathDataInterimTablesCleanedOrganismTranslatedInterim")
    val outputWriter = TsvWriter(GZIPWrite(pathDataInterimTablesCleanedOrganismTranslatedInterim), TsvWriterSettings())
    val headers = dataCleanedOriginalOrganism.firstOrNull()?.keys ?: setOf(
        "organismValue", "value_min", "dbQuality", "value_max",
        "rank", "ids", "taxonomy", "organismCleaned", "organismDbTaxo", "taxonId", "organismInterim"
    )
    outputWriter.writeHeaders(headers)
    dataCleanedOriginalOrganism.forEach { row ->
        if (!row["organismInterim"].isNullOrEmpty())
            outputWriter.writeRow(headers.map { row.getOrDefault(it, "") })
    }
    outputWriter.close()

    logger.info("Writing GNFinder files in $pathDataInterimTablesTranslatedOrganism")
    dataCleanedOriginalOrganism.filter { !it["organismInterim"].isNullOrEmpty() }
        .chunked(cut).mapIndexed { idx, chunk ->
            val number = (idx * cut + cut).toString().padStart(6, '0')
            val fileName = "$pathDataInterimTablesTranslatedOrganism/$number.tsv"
            val writer = TsvWriter(File(fileName), TsvWriterSettings())
            writer.writeHeaders("organismInterim")
            chunk.map { writer.writeRow(it) }
            writer.close()
        }

    logger.info("Completely done. Go Kotlin \\o/")
}