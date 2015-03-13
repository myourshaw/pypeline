__author__ = 'myourshaw'

import mysql.connector

config = {
    'host': 'cortex.local',
    'database': 'vax_test',
    'user': 'vax',
    'password': 'vax',
}

cnx = mysql.connector.connect(**config)
cnx.close()

