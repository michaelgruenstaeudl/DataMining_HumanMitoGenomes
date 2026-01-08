process bwa_index {
    label "bwa_env"

    tag "bwa-index"

    input:
    path ref

    output:
    path "${ref}.*"

    script:
    """
    bwa-mem2 index ${ref}
    """
}
