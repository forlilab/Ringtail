from mpmanager import MPManager


class ResultsManager():

    def __init__(self,
                 mode='dlg',
                 dbman=None,
                 filelist=None,
                 chunk_size=1000,
                 numclusters=3,
                 no_print_flag=False):
        self.dbman = dbman
        self.filelist = filelist
        self.num_result_files = len(filelist)
        self.no_print_flag = no_print_flag
        if mode == 'dlg':
            self.parser = MPManager(filelist=self.filelist,
                                    db_obj=self.dbman,
                                    chunksize=chunk_size,
                                    mode='dlg',
                                    numclusters=numclusters,
                                    no_print_flag=self.no_print_flag)
        # Future parser for vina
        # elif mode == 'vina':
        #     self.parser = MPManager(filelist=self.filelist,
        #                             db_obj=self.dbman,
        #                             chunksize=chunk_size,
        #                             mode='vina',
        #                             numclusters=numclusters,
        #                             no_print_flag=self.no_print_flag)

    def process_results(self):
        # start MP process
        self.parser.process_files()
