#!/usr/bin/env python
# -*- coding: utf-8 -*-

import mysql.connector
from mysql.connector import errorcode
from pprint import pprint

#http://geert.vanderkelen.org/connectorpython-custom-cursors/
class MySQLCursorDict(mysql.connector.cursor.MySQLCursor):
    def _row_to_python(self, rowdata, desc=None):
        row = super(MySQLCursorDict, self)._row_to_python(rowdata, desc)
        if row:
            return dict(zip(self.column_names, row))
        return None

#http://ianhowson.com/a-quick-guide-to-using-mysql-in-python.html
#applies only to MySQLdb?
#see also cursorclass=MySQLdb.cursors.DictCursor
#and http://eric.lubow.org/2010/python/when-to-use-mysql-cursor-classes-in-python/
#for use of built-in cursorclass=MySQLdb.cursors.DictCursor
def FetchOneAssoc(cursor) :
    data = cursor.fetchone()
    if data == None :
        return None
    desc = cursor.description
    dict = {}
    for (name, value) in zip(desc, data) :
        dict[name[0]] = value
    return dict

config = {
  'user': 'test',
  'password': '',
  'host': 'cortex.local',
  'database': 'test',
  'get_warnings': True,
}

#cursor
try:
    cnx = mysql.connector.connect(**config)
    #cnx = mysql.connector.connect(
    #    user='test',
    #    password='',
    #    host='cortex.local',
    #    database='test')
except mysql.connector.Error as err:
  if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
    print("Something is wrong with your user name or password")
  elif err.errno == errorcode.ER_BAD_DB_ERROR:
    print("Database does not exists")
  else:
    print(err)
else:
    cursor = cnx.cursor()
    query = ("SELECT name, address, ssn FROM bar "
             "WHERE name = '%'")
    this_name = 'foobar'
    cursor.execute(query, (this_name))
    for (name, address, ssn) in cursor:
        print '{}\t{}\t{}'.format(name, address, ssn)
finally:
    cursor.close
    cnx.close()

#http://geert.vanderkelen.org/connectorpython-custom-cursors/
#Fetching rows as dictionaries with MySQL Connector
try:
    cnx = mysql.connector.connect(**config)
except mysql.connector.Error as err:
  if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
    print("Something is wrong with your user name or password")
  elif err.errno == errorcode.ER_BAD_DB_ERROR:
    print("Database does not exists")
  else:
    print(err)
else:
    cur = cnx.cursor(cursor_class=MySQLCursorDict)
    cur.execute("SELECT c1, c2 FROM t1")
    rows = cur.fetchall()
    pprint(rows)
    cur.close()
    cnx.close()

finally:
    cursor.close
    cnx.close()
