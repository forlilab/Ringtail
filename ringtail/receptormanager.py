

class ReceptorManager():

    def __init__(self, rec_file_list, dbmanager):

        self.file_list = rec_file_list
        self.receptors = []
        self.dbman = dbmanager

        self._make_receptor_objects()

    def _make_receptor_objects(self):
        for rec_file in self.file_list:
            rec_name = rec_file.split(".")[0]
            rec_obj = make_receptor_object(rec_file, rec_name)
            self.receptors.append(rec_obj)

    def add_receptors_to_db(self):
        for rec in self.receptors:
            rec_name = rec.name
            self.dbman.add_receptor_object_to_row(rec, rec_name)
