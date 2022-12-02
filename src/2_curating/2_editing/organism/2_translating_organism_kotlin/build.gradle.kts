import org.jetbrains.kotlin.gradle.tasks.KotlinCompile

val univocityParserVersion = "2.9.1"
val junitApiVersion = "5.9.1"

group = "net.nprod.onpdb"
version = "0.5-SNAPSHOT"

plugins {
    kotlin("jvm") version "1.7.22"
    application
    id("com.github.johnrengelman.shadow") version "7.1.2"
}

repositories {
    mavenCentral()
    maven {
        url = uri("https://dl.bintray.com/kotlin/kotlin-eap")
    }
}

dependencies {
    implementation(kotlin("stdlib-jdk8"))

    // Univocity is a fast TSV reader/writer
    implementation("com.univocity", "univocity-parsers", univocityParserVersion)
    implementation("org.junit.jupiter", "junit-jupiter", junitApiVersion)

    testImplementation("org.junit.jupiter", "junit-jupiter-api", junitApiVersion)
    testImplementation("org.junit.jupiter", "junit-jupiter", junitApiVersion)
}

/**
 * We need to compile for the outdated conda java version
 */
tasks.withType<KotlinCompile>() {
    kotlinOptions.jvmTarget = "16"
}

/**
 * We say what needs to be executed
 */
application {
    mainClass.set("MainKt")
}

tasks.withType<Jar> {
    manifest {
        attributes(
            mapOf(
                "Main-Class" to application.mainClass
            )
        )
    }
}

/**
 *  We want to make a small jar, with always the same name
 */
tasks.withType<com.github.jengelman.gradle.plugins.shadow.tasks.ShadowJar> {
    minimize()
    archiveBaseName.set("shadow")
    archiveClassifier.set("")
    archiveVersion.set("")

}

/**
 * Configuration of test framework
 */

tasks.test {
    useJUnitPlatform()
}
