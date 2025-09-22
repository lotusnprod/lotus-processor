import org.jetbrains.kotlin.gradle.tasks.KotlinCompile
import org.jetbrains.kotlin.gradle.dsl.JvmTarget

val univocityParserVersion = "2.9.1"
val junitApiVersion = "5.13.4"
val jvmVersion = JvmTarget.JVM_17

group = "net.nprod.onpdb"
version = "0.5-SNAPSHOT"

plugins {
    kotlin("jvm") version "2.2.20"
    application
    id("com.gradleup.shadow") version "9.1.0"
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
tasks.shadowJar {
    minimize()
    archiveBaseName.set("shadow")
    archiveClassifier.set("")
    archiveVersion.set("")
}

/** Configure test framework */
tasks.test {
    useJUnitPlatform()
}