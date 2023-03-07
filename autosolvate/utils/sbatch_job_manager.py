import os
import sys
import time
import subprocess
import logging

from logging import DEBUG, INFO, WARN, WARNING, CRITICAL
logging.basicConfig(level = INFO, force = True, handlers=[])

logger              = logging.getLogger(name = "SbatchJob")
# output_handler      = logging.FileHandler(filename = "log.txt", mode = "a", encoding="utf-8")
output_handler      = logging.StreamHandler()
output_formater     = logging.Formatter(fmt = '%(asctime)s %(name)s %(levelname)s: %(message)s', datefmt="%H:%M:%S")
output_handler.setFormatter(output_formater)
logger.addHandler(output_handler)


class SbatchJob():
    def __init__(self, script:str, maxtime = 86400*7, submit_dir = "", check_interaval = 10):
        """A data class that stores a slurm job information"""
        if not isinstance(script, str):
            raise ValueError("script file path required. the input is a "+ str(type(script)))
        if not os.path.exists(script):
            raise OSError("file not found: " + script)
        self.name = script
        with open(script, "r") as f:
            self.script = f.read()
        self.jobid          = -1
        self.partition      = ""
        self.time           = 0
        self.jobname        = ""
        self.user           = ""
        self.status         = "U"
        self.nodelist       = ""
        self.output         = ""
        if not submit_dir:
            submit_dir = os.path.dirname(script)
            if not submit_dir:
                submit_dir = os.getcwd()
        self.submit_dir = submit_dir

        self.maximum_time = maxtime
        self.check_interval = check_interaval



    def asctime2second(self, s:str):
        if s.count(":") == 1:
            sec, miu = int(s[-2:]), int(s[:-3])
            self.time = sec + 60 * miu
        elif s.count(":") == 2 and s.count("-") == 0:
            sec, miu, hor = int(s[-2:]), int(s[-5:-3]), int(s[:-6])
            self.time = sec + 60 * miu + 3600 * hor
        elif s.count(":") == 2 and s.count("-") == 1:
            d, s = s.split("-")
            sec, miu, hor = int(s[-2:]), int(s[-5:-3]), int(s[:-6])
            self.time = sec + 60 * miu + 3600 * hor + 86400 * int(d)
        else:
            logger.warn("warning: the time " + s + " cannot be processed.")
            self.time = self.time
        return self.time

    def check_finished(self):
        data = subprocess.getoutput("squeue").splitlines()
        if len(data) == 1:
            self.status = "F"
            return 
        for line in data[1:]:
            args = line.split()
            jid = int(args[0])
            if jid == self.jobid:
                self.partition  = args[1]
                self.jobname    = args[2]
                self.user       = args[3]
                self.status     = args[4]
                self.time       = self.asctime2second(args[5])
                self.nodelist   = args[-1]
                break
        else:
            self.status = "F"

    def submit(self):
        cwd = os.getcwd()
        if not os.path.exists(self.submit_dir):
            os.makedirs(self.submit_dir)
        os.chdir(self.submit_dir)
        result = subprocess.getoutput("sbatch " + self.name)
        os.chdir(cwd)
        jobid = result.split()[-1]
        try:
            self.jobid = int(jobid)
            self.output = os.path.join(os.getcwd(), "slurm-"+jobid+".out")
            self.check_finished()
        except Exception:
            logger.warn("failed to submit script " + self.name + ". slurm report:")
            logger.warn(result)

    def wait_until_finish(self):
        if self.status == "U":
            self.submit()
        for i in range(self.maximum_time // self.check_interval):
            time.sleep(self.check_interval)
            self.check_finished()
            if self.status == "F":
                logger.info(f"job {self.jobid:9d} {self.jobname:9s} finished. ")
                break
            else:
                logger.info(f"job {self.jobid:9d} {self.jobname:9s} status {self.status:2s} time {self.time:9d}")
        else:
            logger.warn(f"job {self.jobid:9d} {self.jobname:9s} has been running for more than {self.maximum_time} seconds, no longer waiting for this task to complete.")


class SbatchJobManager():
    def __init__(self, scriptlist:list = [], maxtime = 86400*7, check_interval = 10):
        if isinstance(scriptlist, str):
            scriptlist = [scriptlist, ]
        self.scriptlist = scriptlist
        self.joblist = []
        self.jobstatus = []
        for script in scriptlist:
            self.joblist.append(SbatchJob(script, maxtime = maxtime))
            self.jobstatus.append(self.joblist[-1].status)
        self.maxtime = maxtime
        self.check_interval = check_interval
        self.allfinished = False

    def append(self, sjob:SbatchJob):
        if isinstance(sjob, str):
            sjob = SbatchJob(sjob, maxtime=self.maxtime)
        self.joblist.append(sjob)
        self.jobstatus.append(sjob.status)
        
    def submit_all(self):
        for sjob in self.joblist:
            sjob:SbatchJob
            if sjob.status == "U":
                sjob.submit()
                print(f"job {sjob.jobid:9d} {sjob.jobname:9s} submitted. ")

    def wait_until_finish(self):
        self.submit_all()
        for i in range(self.maxtime // self.check_interval):
            time.sleep(self.check_interval)
            n_running = 0
            for i, sjob in enumerate(self.joblist):
                if sjob.status == "F" and self.jobstatus[i] == "F":
                    continue
                sjob.check_finished()
                self.jobstatus[i] = sjob.status
                if sjob.status == "F":
                    print(f"job {sjob.jobid:9d} {sjob.jobname:9s} finished. ")
                else:
                    print(f"job {sjob.jobid:9d} {sjob.jobname:9s} status {sjob.status:2s} time {sjob.time:9d}")
                    n_running += 1
            if n_running == 0:
                print(f"all {len(self.joblist)} jobs are finished. total time {i * self.check_interval} seconds.")
                break
            else:
                print(f"{n_running} out of {len(self.joblist)} jobs are still running.")
        else:
            n_running = 0
            for i, sjob in enumerate(self.joblist):
                if sjob.status != "F":
                    n_running += 1
                    print(f"job {sjob.jobid:9d} {sjob.jobname:9s} has been running for more than {self.maxtime} seconds, no longer waiting for this task to complete.")
        print(f"{len(self.joblist)-n_running} out of {len(self.joblist)} jobs have finished.")



