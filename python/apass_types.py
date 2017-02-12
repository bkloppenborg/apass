import sys
import os

# Custom modules
sys.path.append(os.path.abspath('./modules/FileLock/filelock'))
from filelock import FileLock

class FileStore():

    def __init__(self, filename):
        self.filename = filename

    def close(self):
        self.file.close()
        self.lock.release()

    def open(self):
        # obtain a lock or keep retrying for 60 seconds
        self.lock = FileLock(self.filename, timeout=60, delay=0.05)
        self.lock.acquire()
        self.file = open(self.filename, 'a+')

    def append(self, data):
        self.file.write(data)

class FileStores():
    files = {}

    def __init__(self, datapath):
        self.datapath = datapath

    def close(self):
        for id, filestore in self.files.iteritems():
            filestore.close()

    def get_filehandle(self, file_id):

        if file_id not in self.files:
            filename = self.datapath + "/" + str(file_id).zfill(5) + ".dat"
            filestore = FileStore(filename)
            filestore.open()
            self.files[file_id] = filestore

        return self.files[file_id]

    def insert(self, file_id, data):
        filehandle = self.get_filehandle(file_id)
        filehandle.append(data)
