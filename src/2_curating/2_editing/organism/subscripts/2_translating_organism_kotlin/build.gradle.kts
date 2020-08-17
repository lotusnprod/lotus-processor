import org.jetbrains.kotlin.gradle.tasks.KotlinCompile

val univocityParserVersion = "2.8.4"
group = "net.nprod.onpdb"
version = "0.3-SNAPSHOT"

plugins {
    kotlin("jvm") version "1.3.72"
    application
    id("com.github.johnrengelman.shadow") version "5.0.0"
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
}

/**
 * We need to compile for the outdated conda java version
 */
tasks.withType<KotlinCompile>() {
    kotlinOptions.jvmTarget = "1.8"
}

/**
 * We say what needs to be executed
 */
application {
    mainClassName = "MainKt"
}

tasks.withType<Jar> {
    manifest {
        attributes(
            mapOf(
                "Main-Class" to application.mainClassName
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