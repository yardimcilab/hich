onsuccess:
    """echo {config[email_body_success]:q} | mail -s {config[email_subject_success]:q} {config[email_address]:q}"""

onerror:
    """echo {config[email_body_error]:q} | mail -s {config[email_subject_error]:q} {config[email_address]:q}"""