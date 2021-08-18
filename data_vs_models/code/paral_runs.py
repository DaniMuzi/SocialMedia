#!/usr/bin/python

import sys
import time
import subprocess
import psutil


def pid_is_running(pid):
    p = psutil.Process(pid)
    if p.status() == psutil.STATUS_ZOMBIE:
        return False
    return True



dataset = sys.argv[1]
delta = sys.argv[2]
Smin = sys.argv[3]

NumFiles = int(sys.argv[4])
NumCores = int(sys.argv[5])



PIDs = []

for num in range(NumFiles):
  process = subprocess.Popen(['./data_vs_models.out', dataset, delta, Smin, str(num+1)])
  PIDs.append(process.pid)

  while len(PIDs) == NumCores:

      time.sleep(60)

      done = []
      for i, p in enumerate(PIDs):
          if pid_is_running(p) == False:
              done.append(i)

      for i in reversed(done):
          PIDs.pop(i)


