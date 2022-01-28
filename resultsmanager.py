from mpmanager import MPManager

class ResultsManager():

    def __init__(self, mode='dlg', dbman=None, filelist=None, chunk_size=1000, numclusters = 3):
        self.dbman = dbman
        self.filelist = filelist
        if mode =='dlg':
            self.parser = MPManager(filelist = self.filelist, db_obj = self.dbman, chunksize = chunk_size, mode='dlg', numclusters = numclusters)
        #elif mode == 'vina'
            #self.parser = fileparsers.VinaParser()

    def process_results(self):
        #start MP process
        self.parser.process_files()