import subprocess
import time

def submit_job(command):
    print("[run command]: ", command)
    while True:
        try:
            output = subprocess.check_output(command, shell=True, text=True)
            print(output)
            break
        except subprocess.CalledProcessError as e:
            print(f"error: {e}")
            time.sleep(1)
    time.sleep(0.5)

def submit_job_once(command):
    print("[run command]: ", command)
    try:
        output = subprocess.check_output(command, shell=True, text=True)
        print(output)
    except subprocess.CalledProcessError as e:
        print(f"error: {e}")
        time.sleep(1)
