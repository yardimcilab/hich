onstart:
    shell("""echo \"Success/error email alert with subject [Snakemake Alert] will be sent to {config[email_address]}\"""")

onsuccess:
    shell("""echo \"{config[email_body_success]}\" | mail -s \"{config[email_subject_success]}\" \"{config[email_address]}\"""")

onerror:
    shell("""echo \"{config[email_body_error]}\" | mail -s \"{config[email_subject_error]}\" \"{config[email_address]}\"""")
