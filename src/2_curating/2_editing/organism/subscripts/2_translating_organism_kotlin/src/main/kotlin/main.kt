import com.univocity.parsers.tsv.TsvWriter
import com.univocity.parsers.tsv.TsvWriterSettings

import java.io.File
import java.text.DecimalFormat
import java.util.stream.Collectors
import kotlin.system.exitProcess

// By how much the gnfinder files are going to be cut into

const val cut = 10_000

// Show usage and exit
fun usageExit() {
    println("Usage: java -jar <scriptname> Data_path full/min")
    println("If the last argument is min, it will run on tables if not on tables_min")
    println("If the last argument is full, it will run on tables")
    exitProcess(1)
}

fun main(args: Array<String>) {
    val logger = MyDirtyLogger() //LoggerFactory.getLogger("main")
    if (args.size < 2) usageExit()
    if (args[1] !in setOf("min", "full")) usageExit()

    val pathData = args[0]

    val pathDataInterim = "$pathData/interim"

    val pathDataInterimDictionaries = "$pathDataInterim/dictionaries"
    val pathDataInterimDictionariesTcmNames = "$pathDataInterimDictionaries/tcm/names.tsv.gz"
    val pathDataInterimDictionariesCommonNames = "$pathDataInterimDictionaries/common/names.tsv.gz"
    val pathDataInterimDictionariesCommonBlackDic = "$pathDataInterimDictionaries/common/black.tsv"

    val pathDataInterimTables = when (args[1]) {
        "full" -> "$pathData/interim/tables"
        "min" -> "$pathData/interim/tables_min"
        else -> throw Exception("This shouldn't have happened, we only know these types of tables")
    }

    val pathDataInterimTablesCleaned = "$pathDataInterimTables/2_cleaned"
    val pathDataInterimTablesCleanedOrganism = "$pathDataInterimTablesCleaned/organism"
    val pathDataInterimTablesCleanedOrganismOriginalTable = "$pathDataInterimTablesCleanedOrganism/original.tsv.gz"
    val pathDataInterimTablesCleanedOrganismTranslatedInterim =
        "$pathDataInterimTablesCleanedOrganism/interim.tsv.gz"

    val pathDataInterimTablesTranslatedOrganism = "$pathDataInterimTables/1_translated/organism"

    logger.info("Making sure directories exist")
    File(pathDataInterimTablesTranslatedOrganism).mkdirs()
    File(pathDataInterimTablesCleanedOrganism).mkdirs()

    // Organism List

    logger.info("Loading and processing the organism list")
    val dataCleanedOriginalOrganism = parseTSVZFile(pathDataInterimTablesCleanedOrganismOriginalTable)?.map {
        it.toFieldMap()
    } ?: throw Exception("Sorry can't read organism list.")

    // TCM

    logger.info("Loading and processing the TCM names")

    val tcmNamesDic = parseTSVZFile(pathDataInterimDictionariesTcmNames)?.flatMap {
        val vernacularName = it.getValue("vernacularName", "")
        val canonicalName = it.getValue("canonicalName", "")
        val newCanonicalName = it.getValue("newCanonicalName", "")
        val out = mutableListOf(
            Pair(Regex("\\b${vernacularName}\\b", RegexOption.IGNORE_CASE), canonicalName)
        )
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
    val exclusionDic = parseTSVFile(pathDataInterimDictionariesCommonBlackDic)?.map {
        if (it.getString(0) == "Apple") println("Apple will become empty")
        it.getValue<String>("blackName", null)
    }?.sortedByDescending { it.length }?.map { Regex("\\b$it\\b", RegexOption.IGNORE_CASE) }?.distinct()
        ?: throw Exception("Sorry can't read exclusion list")

    // Processing

    logger.info("Processing")
    val processedRecords =
        dataCleanedOriginalOrganism.parallelStream().map { record -> // We remove the ones from the exclusion list

            var newTerm = record["organismOriginal"] ?: ""
            exclusionDic.map {
                newTerm = newTerm.replace(it, "")
            }

            newTerm = newTerm.replace(
                Regex("\\b${record.getOrDefault("organismCleaned", "")}\\b", RegexOption.IGNORE_CASE),
                ""
            )

            // We try to clean up a bit

            newTerm = newTerm.replace(".", "").trim().replace(" ", " ").replace("_", " ")
                .replace(Regex(" +"), " ")

            commonNamesDic.forEach {
                newTerm = newTerm.replace(it.first, it.second)
            }

            tcmNamesDic.forEach {
                newTerm = newTerm.replace(it.first, it.second)
            }

            // Second round of TCM

            tcmNamesDic.forEach {
                newTerm = newTerm.replace(it.first, it.second)
            }

            record["organismInterim"] = newTerm
            record
        }.collect(Collectors.toList())

    logger.info("Done processing")

    // Writing files

    logger.info("Writing Interim file $pathDataInterimTablesCleanedOrganismTranslatedInterim")
    val outputWriter = TsvWriter(GZIPWrite(pathDataInterimTablesCleanedOrganismTranslatedInterim), TsvWriterSettings())
    val headers = processedRecords.first()?.keys ?: throw Exception("We don't even have a header to write.")
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