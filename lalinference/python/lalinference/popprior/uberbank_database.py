'''
For N templates in the uberbank, the overlaps form an NxN symmetric matrix. Currently, there are many subbank overlap files that comprise the uberbank, each of which is a rectangular section of the non-symmetric part of the overlap matrix. The subbanks are essentially arbitrary in size. Templates with overlaps < 0.25 are also not computed (this may need to change in the future).
This code creates a database for locating the overlaps between templates in the uberbank. This database allows one to retrieve all the overlaps for a particular template. The overlaps are placed into a 1xN numpy array and are ordered.
'''

import glob
import h5py
import sqlite3
import time
from glue.text_progress_bar import ProgressBar
import itertools
import numpy as np
import os

class Bank(object):
    def __init__(self, connection):
        self.connection = connection
        cursor = self.connection.cursor()
        self.numtemplates, = cursor.execute("SELECT MAX(k)+1 FROM bank;").fetchall() # finds total number of templates in bank

    @staticmethod
    def make_db(db_filename, filenames):
        # Make the sqlite database if one does not exist
        connection = sqlite3.connect(db_filename)
        cursor = connection.cursor()
        # Create Table 1 that contains:
        # - the filenames of the hdf5 subbank files
        # - their respective ID integers
        cursor.execute("CREATE TABLE fname (filename_id INTEGER PRIMARY KEY, filename TEXT);")
        # Create Table 2 that contains:
        # - k (integers which will be populated later, this column gives the order in which overlaps are placed),
        # - ID integers of the filenames (instead of saving the entire filename to the database - this is to save memory)
        # - m1, m2, chi1, chi2 of templates (separate column for each)
        # - row and column in which this overlap is found in the file (if row=None, means that the overlaps are found in the columns only)
        cursor.execute("CREATE TABLE bank (k INTEGER, filename_id INTEGER, m1 REAL, m2 REAL, chi1 REAL, chi2 REAL, row INTEGER, column INTEGER);")
        # Populate Table 1:
        for filename in filenames:
            cursor.execute("INSERT INTO fname (filename) VALUES (?)", (filename,))
        connection.commit()
        cursor.execute("CREATE INDEX filename_index ON fname (filename)")
        progress = ProgressBar(max=len(filenames))
        # Populate Table 2:
        for filename in filenames:
            progress.increment(text=filename)
            try:
                f = h5py.File(filename,"r")
                nrows = f[f.keys()[0]]['overlaps'].shape[0] # find number of rows in subbank file
                for i, (m1, m2, chi1, chi2) in enumerate(zip(f[f.keys()[0]]['mass1'].value, f[f.keys()[0]]['mass2'].value, f[f.keys()[0]]['spin1z'].value, f[f.keys()[0]]['spin2z'].value)):
                    cursor.execute("INSERT INTO bank (filename_id, m1, m2, chi1, chi2, row, column) VALUES ((SELECT filename_id FROM fname WHERE filename = ?),?,?,?,?,?,?)", (filename, m1, m2, chi1, chi2, i if i <= nrows else None, i))
                    #FIXME: After building the database, confirm that each template only appears once in one row of the hdf5 files
            except: # ignore corrupted/broken subbank files
                print "Cannot load h5py file:", filename
        cursor.execute("CREATE INDEX template_index ON bank (m1, m2, chi1, chi2)")
        cursor.execute("CREATE INDEX filename_id_index ON bank (filename_id)")
        cursor.execute("CREATE TEMPORARY TABLE template AS SELECT DISTINCT m1,m2,chi1,chi2 FROM bank;")
        cursor.execute("CREATE INDEX tmp ON template (m1,m2,chi1,chi2);")
        cursor.execute("UPDATE bank SET k = (SELECT rowid-(SELECT MIN(rowid) FROM template) FROM template WHERE m1 = bank.m1 and m2 = bank.m2 and chi1 = bank.chi1 and chi2 = bank.chi2);") # populate column k
        cursor.execute("DROP TABLE template;")
        cursor.execute("CREATE INDEX k_index ON bank (k);")
        connection.commit()
        return connection

    def get_overlaps(self, template_number):
        # Obtain the overlaps, given a template_number (integer)
        overlaps = np.zeros(self.numtemplates)
        cursor = self.connection.cursor()
        for (filename, column), k in itertools.groupby(cursor.execute("SELECT fname.filename, a.column, b.k FROM bank AS a JOIN bank AS b ON (b.filename_id = a.filename_id) JOIN fname ON (fname.filename_id=a.filename_id) WHERE a.k = ? AND a.row IS NULL AND b.row IS NOT NULL ORDER BY fname.filename, b.row;", (template_number,)), lambda (fn, c, k_idx): (fn, c)):
            f = h5py.File(filename, "r")
            overlaps[[x[-1] for x in k]] = f[f.keys()[0]]['overlaps'].value[:,column]
        for (filename, row), k in itertools.groupby(cursor.execute("SELECT fname.filename, a.row, b.k FROM bank AS a JOIN bank AS b ON (b.filename_id = a.filename_id) JOIN fname ON (fname.filename_id = a.filename_id) WHERE a.k = ? AND a.row IS NOT NULL ORDER BY fname.filename, b.column;", (template_number,)), lambda (fn, r, k_idx): (fn, r)):
            f = h5py.File(filename, "r")
            overlaps[[x[-1] for x in k]] = f[f.keys()[0]]['overlaps'].value[row,:]
        return overlaps

    def get_templates(self):
        # Obtain the parameters (m1, m2, chi1, chi2) of the templates, in order k
        cursor = self.connection.cursor()
        return cursor.execute("SELECT DISTINCT m1,m2,chi1,chi2 FROM bank ORDER BY k;").fetchall()


'''The following files currently cannot be read'''
#Cannot load h5py file: uberbank/bank_299_overlaps.hdf
#Cannot load h5py file: uberbank/bank_479_overlaps.hdf
#Cannot load h5py file: uberbank/bank_527_overlaps.hdf
#Cannot load h5py file: uberbank/bank_528_overlaps.hdf
#Cannot load h5py file: uberbank/bank_529_overlaps.hdf
#Cannot load h5py file: uberbank/bank_548_overlaps.hdf
#Cannot load h5py file: uberbank/bank_549_overlaps.hdf
#Cannot load h5py file: uberbank/bank_559_overlaps.hdf
#Cannot load h5py file: uberbank/bank_579_overlaps.hdf
#Cannot load h5py file: uberbank/bank_589_overlaps.hdf
#Cannot load h5py file: uberbank/bank_619_overlaps.hdf
#Cannot load h5py file: uberbank/bank_688_overlaps.hdf
#Cannot load h5py file: uberbank/bank_689_overlaps.hdf
#Cannot load h5py file: uberbank/bank_739_overlaps.hdf
#Cannot load h5py file: uberbank/bank_839_overlaps.hdf
        
start_time = time.time()
#Bank.make_db("uberbank_database.sqlite",glob.glob("uberbank/*hdf")) # make sqlite file
x = Bank(sqlite3.connect("uberbank_database.sqlite"))
#x = Bank(Bank.make_db(":memory:", glob.glob("uberbank/*.hdf")) # in RAM version
print "Seconds taken to create database:", time.time()-start_time
