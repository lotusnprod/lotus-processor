tasks.register("castShadows") {
    dependsOn(gradle.includedBuild("2_translating_organism_kotlin").task(":shadowJar"))
}
