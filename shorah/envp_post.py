import gzip

def post_process_envp(template_file, file_to_process):
    with gzip.open(template_file, "rt") as template, gzip.open(file_to_process, "wt") as f:
        pass