
import inspect
import multiprocessing
import sys

import pyRootPwa.utils

stdoutLock = multiprocessing.Lock()
stdoutisatty = None
stderrisatty = None

__terminalColorStrings = {}
__terminalColorStrings['normal']     = "\033[0m"
__terminalColorStrings['bold']       = "\033[1m"
__terminalColorStrings['underline']  = "\033[4m"
__terminalColorStrings['blink']      = "\033[5m"
__terminalColorStrings['inverse']    = "\033[7m"
__terminalColorStrings['fgBlack']    = "\033[30m"
__terminalColorStrings['fgRed']      = "\033[31m"
__terminalColorStrings['fgGreen']    = "\033[32m"
__terminalColorStrings['fgYellow']   = "\033[33m"
__terminalColorStrings['fgBlue']     = "\033[34m"
__terminalColorStrings['fgMangenta'] = "\033[35m"
__terminalColorStrings['fgCyan']     = "\033[36m"
__terminalColorStrings['fgWhite']    = "\033[37m"
__terminalColorStrings['bgBlack']    = "\033[40m"
__terminalColorStrings['bgRed']      = "\033[41m"
__terminalColorStrings['bgGreen']    = "\033[42m"
__terminalColorStrings['bgYellow']   = "\033[43m"
__terminalColorStrings['bgBlue']     = "\033[44m"
__terminalColorStrings['bgMangenta'] = "\033[45m"
__terminalColorStrings['bgCyan']     = "\033[46m"
__terminalColorStrings['bgWhite']    = "\033[47m"

def __printFormatted(msg, level):
	frame = inspect.currentframe().f_back.f_back
	if frame is None:
		printErr("This method cannot be called directly.")
	(filename, lineno, function, code_contex, index) = inspect.getframeinfo(frame)
	if function == "<module>":
		function = "__main__"
	string = ""
	if level == "err":
		if pyRootPwa.utils.stderrisatty: string += __terminalColorStrings['fgRed']
		string += "!!! "
		string += function + " [" + filename + ":" + str(lineno) + "]: "
		string += "error: "
		if pyRootPwa.utils.stderrisatty: string += __terminalColorStrings['normal']
		string += msg
	elif level == "warn":
		if pyRootPwa.utils.stderrisatty: string += __terminalColorStrings['fgYellow']
		string += "??? "
		string += function + " [" + filename + ":" + str(lineno) + "]: "
		string += "warning: "
		if pyRootPwa.utils.stderrisatty: string += __terminalColorStrings['normal']
		string += msg
	elif level == "suc":
		if pyRootPwa.utils.stdoutisatty: string += __terminalColorStrings['fgGreen']
		string += "*** "
		string += function + ": success: "
		if pyRootPwa.utils.stdoutisatty: string += __terminalColorStrings['normal']
		string += msg
	elif level == "info":
		if pyRootPwa.utils.stdoutisatty: string += __terminalColorStrings['bold']
		string += ">>> "
		string += function + ": info: "
		if pyRootPwa.utils.stdoutisatty: string += __terminalColorStrings['normal']
		string += msg
	elif level == "debug":
		if pyRootPwa.utils.stdoutisatty: string += __terminalColorStrings['fgMangenta']
		string += "+++ "
		string += function + ": debug: "
		string += __terminalColorStrings['normal']
		string += msg
	else:
		printErr("Invalid level string.")
		raise Exception()
	if level == "err" or level == "warn":
		sys.stderr.write(string + "\n")
	else:
		sys.stdout.write(string + "\n")


def printErr(msg): __printFormatted(str(msg), "err")
def printWarn(msg): __printFormatted(str(msg), "warn")
def printSucc(msg): __printFormatted(str(msg), "suc")
def printInfo(msg): __printFormatted(str(msg), "info")
def printDebug(msg): __printFormatted(str(msg), "debug")

