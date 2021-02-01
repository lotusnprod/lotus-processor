fun String.clean() = this.replace(Regex(" +"), " ").trim()

data class Dictionaries(
    val exclusionDic: MutableList<Regex> = mutableListOf(),
    val combinedDic: MutableList<Pair<Regex, String>> = mutableListOf(),
) {
    /**
     * This function, will take a record and return a new organismInterim with the names from the dictionary
     * : removed for exclusionDic and replaced for combinedDic
     */
    @Suppress("NOTHING_TO_INLINE")
    inline fun processRecord(record: Map<String, String>): String {
        var newTerm = record["organismValue"]!!.replace(Regex("[. _]"), " ")

        // We process the term with the exclusion dictionary

        exclusionDic.map {
            newTerm = newTerm.replace(it, "")
        }

        // If we already removed everything we exit

        if (newTerm == "") return ""

        // If we have an organism Cleaned, we check if we are not equal, in which case we have nothing to do
        // if not, we remove it from the string anyway to not have duplicates

        val orgCleaned = record["organismCleaned"]!!
        if (orgCleaned == newTerm) return ""


        newTerm = newTerm.replace(orgCleaned, "", ignoreCase = true)


        // We try to clean up a bit

        newTerm = newTerm.clean()

        if (newTerm == "") return ""

        val candidateList = mutableSetOf<String>()
        // So this is a bit complicated
        // We match each string only once and remove it from the initial string so it doesn't get matched again
        //
        // So for "foo bar bim" if we have in the dictionary "foo bar" => ABC  and "bim" => XYZ and "foo" => 123
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
                    // We remove it, find the index of the matched string as lowercase and remove it
                    newTerm = newTerm.replace(pair.first, "")
                    // And add the replacement to the candidateList
                    if (pair.second != "")
                        candidateList.add(pair.second)
                    if (pair.second == "Pagellus erythrinus")
                        println("What gave it is $pair")
                    newTerm.clean()

                    if (newTerm == "") return@loop
                }
            }
        }

        // We combine both the new term and the eventual candidateList
        return candidateList.joinToString(" , ")
    }
}