from mpmanager import MPManager

class ResultsManager():
    
    def __init__(self, mode='dlg', dbman=None, filelist=None, chunk_size=1000, numclusters):
        self.dbman = dbman
        self.filelist = filelist
        if mode =='dlg':
            self.dlg_parser = MPManager(filelist = self.filelist, db_obj = self.dbman, chunksize = chunk_size, mode='dlg')
        #elif mode == 'vina'
            #self.parser = fileparsers.VinaParser()

        self.process_results()

    def process_results(self):
        #start MP process
        self.dlg_parser.process_files()