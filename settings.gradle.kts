pluginManagement {
    repositories {
        mavenCentral()
        gradlePluginPortal()
        maven {
            url = uri("https://dl.bintray.com/kotlin/kotlin-eap")
        }
    }

}
rootProject.name = "onpdb_scripts"


includeBuild("src/2_curating/2_editing/organism/2_translating_organism_kotlin")