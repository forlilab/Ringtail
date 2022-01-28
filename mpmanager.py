import multiprocessing
from time import sleep
from mpreaderwriter import DockingFileReader
from mpreaderwriter import Writer

class MPManager():

	def __init__(self, filelist, mode, db_obj, chunksize, numclusters):
		self.mode = mode
		self.filelist = filelist
		self.db = db_obj
		self.chunksize = chunksize
		self.numclusters = numclusters

		self.max_proc = multiprocessing.cpu_count()
		self.queueIn = multiprocessing.Queue(maxsize=self.max_proc)
		self.queueOut = multiprocessing.Queue()

	def process_files(self):
	    # start the workers in background
	    for i in range(self.max_proc):
	        # one worker is started for each processor to be used
	        s = DockingFileReader(self.queueIn, self.queueOut, self.mode, self.numclusters)
	        # this method calls .run() internally
	        s.start()
	    # start the writer (only one writer, processing the data from the workers)
	    w = Writer(self.queueOut, self.max_proc, self.chunksize, self.db)
	    w.start()

	    # process items in the queue
	    for file in self.filelist:
	        print("Feeding the queue |%s|" % file)
	        self.queueIn.put(file, block=True)
	    # put as many poison pills in the queue as there are workers
	    for i in range(self.max_proc):
	        self.queueIn.put(None)