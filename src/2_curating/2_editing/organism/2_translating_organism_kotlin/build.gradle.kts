import org.jetbrains.kotlin.gradle.tasks.KotlinCompile
import org.jetbrains.kotlin.gradle.dsl.JvmTarget
import com.github.jengelman.gradle.plugins.shadow.tasks.ShadowJar

val univocityParserVersion = "2.9.1"
val junitApiVersion = "5.13.4"
val jvmVersion = JvmTarget.JVM_17

group = "net.nprod.onpdb"
version = "0.5-SNAPSHOT"

plugins {
    kotlin("jvm") version "2.2.0"
    application
    id("com.github.johnrengelman.shadow") version "8.1.1"
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
tasks.withType<KotlinCompile>().configureEach {
    compilerOptions {
        jvmTarget.set(jvmVersion)
    }
}

tasks.withType<JavaCompile> {
    sourceCompatibility = "17"
    targetCompatibility = "17"
    options.encoding = "UTF-8"
}

/** Application main class */
application {
    mainClass.set("MainKt")
}

/** Add Main-Class to Jar manifest */
tasks.withType<Jar> {
    manifest {
        attributes(
            "Main-Class" to application.mainClass.get()
        )
    }
}

/** ShadowJar configuration */
tasks.withType<ShadowJar> {
    minimize()
    archiveBaseName.set("shadow")
    archiveClassifier.set("")
    archiveVersion.set("")
}

/** Configure test framework */
tasks.test {
    useJUnitPlatform()
}
