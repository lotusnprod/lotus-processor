import com.univocity.parsers.common.record.Record
import com.univocity.parsers.tsv.TsvParser
import com.univocity.parsers.tsv.TsvParserSettings
import java.io.*
import java.time.LocalDateTime
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream

/**
 * Get a BufferedReader from a gzipped file.
 */
fun GZIPRead(name: String): BufferedReader {
    return BufferedReader(InputStreamReader(GZIPInputStream(FileInputStream(name))))
}

/**
 * Get a BufferedReader from a gzipped file.
 */
fun GZIPWrite(name: String): BufferedWriter {
    return BufferedWriter(OutputStreamWriter(GZIPOutputStream(FileOutputStream(name))))
}

/**
 * Get a list of records from the given Reader
 */
fun parseTSVFile(file: Reader): List<Record>? {
    val settingsParser = TsvParserSettings()
    settingsParser.format.setLineSeparator("\n")
    settingsParser.isHeaderExtractionEnabled = true
    val tsvParser = TsvParser(settingsParser)

    return tsvParser.parseAllRecords(file)
}

/**
 * Get a list of records from the filename pointing to a Gzipped TSV file
 */
fun parseTSVZFile(file: String): List<Record>? = parseTSVFile(GZIPRead(file))

/**
 * Get a list of records from the filename pointing to a TSV file
 */
fun parseTSVFile(file: String): List<Record>? = parseTSVFile(File(file).bufferedReader())

/**
 * A really basic logger to avoid having to load slf4j
 */
class MyDirtyLogger {
    /*
        Display `msg` on the console with current time
     */
    fun info(msg: String) = println("${this.time()} - $msg")
    private fun time(): LocalDateTime = LocalDateTime.now()
}
