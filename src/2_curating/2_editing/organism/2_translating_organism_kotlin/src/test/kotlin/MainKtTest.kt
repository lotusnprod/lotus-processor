import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Test

internal class CleaningTest {
    @Test
    fun processRecordOk() {
        val dics = Dictionaries()
        dics.combinedDic.add((Regex("""\bRose Apple\b""", RegexOption.IGNORE_CASE) to "Syzygium jambos"))
        val records = mapOf("organismValue" to "from Rose apple", "organismCleaned" to "")
        val out = dics.processRecord(records)
        assertEquals("Syzygium jambos", out)
    }

    @Test
    fun processRecordDouble() {
        val dics = Dictionaries()
        dics.combinedDic.add((Regex("""\bRose Apple\b""", RegexOption.IGNORE_CASE) to "Syzygium jambos"))
        val out = dics.processRecord(
            mapOf(
                "organismValue" to "from Syzygium, Rose apple",
                "organismCleaned" to "Syzygium"
            )
        )
        assertEquals("Syzygium jambos", out)
    }
}
