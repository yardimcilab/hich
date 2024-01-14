onstart:
    shell("""
    mkdir -p logs/email;
    echo \"Success/error email alert with subject [Snakemake Alert] will be sent to {config[email_address]}\" 2>> logs/email/email.log
    """)

onsuccess:
    shell(
    """
    python3 benchreport/benchreport.py benchmarks benchreport;
    if [ -n "{config[email_address]}" ]; then
        echo "{config[email_body_success]}" | mail -s "{config[email_subject_success]}" "{config[email_address]}" 2>> logs/email/email.log;
    fi
    """)

onerror:
    shell(
    """
    if [ -n "{config[email_address]}" ]; then
        echo "{config[email_body_error]}" | mail -s "{config[email_subject_error]}" "{config[email_address]}" 2>> logs/email/email.log;
    fi
    """)
