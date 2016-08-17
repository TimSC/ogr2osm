# -*- coding: utf-8 -*-

class Feature(object):
	def __init__(self):
		pass

class Point(object):
	def __init__(self, x, y):
		self.x = x
		self.y = y

class Way(object):
	def __init__(self):
		self.points = []

class Relation(object):
	def __init__(self):
		self.members = []

