from dbmanager import DBManager
from resultsmanager import ResultsManager

class VSManager():
    """ DOCUMENTATION GOES HERE """
    def __init__(self, db_fname, db_opts, rman_opts):
        self.dbman = DBManager(db_fname, opts=db_opts)
        self.results_man = ResultsManager(mode=rman_opts['mode'], dbman = self.dbman, chunk_size=rman_opts['chunk_size'], file_list=rman_opts['filelist'], numclusters=rman_opts['num_clusters'])


    def add_results(self):
        """"""
        pass


    def filter(self):
        """"""
        pass




if __name__ == '__main__':
    VManager()
    pass

