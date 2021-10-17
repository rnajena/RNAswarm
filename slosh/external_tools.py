import subprocess

def run_r_script(script_path, *arguments):
    """
    """
    assert script_path[-2:] == ".r", f"{script_path} is not an R script"
    argument_list = ["Rscript", script_path] + arguments
    subprocess.run(argument_list)

def run_shell_script(script_path, *arguments):
    """
    """
    assert script_path[-3:] == ".sh", f"{script_path} is not a shell script"
    argument_list = ["bash", script_path] + arguments
    subprocess.run(argument_list)