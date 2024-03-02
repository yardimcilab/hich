process CAT {
    input: path(csv)
    output: stdout
    script: """cat ${csv}"""
}