#!/usr/bin/env python
# -*- coding: utf-8 -*-

import geom

class MemoryGeom(object):
	def __init__(self):
		self.features=[]
		self.elementIdCounter = -1
		self.elementIdCounterIncr = -1
		self.nodePosDict = {}

	def SetElementIdCounter(self, idnum):
		self.elementIdCounter = idnum
		
	def GetElementIdCounter(self):
		return self.elementIdCounter

	def SetElementIdCounterIncr(self, inc):
		self.elementIdCounterIncr = inc

	def DetermineNodeId(self, geometry, options):
		rx = int(round(geometry.x * 10**(options.significantDigits-options.roundingDigits)))
		ry = int(round(geometry.y * 10**(options.significantDigits-options.roundingDigits)))
		pos = (rx, ry)
		try:
			return self.nodePosDict[pos]
		except KeyError:
			nid = self.elementIdCounter
			self.nodePosDict[pos] = nid
			self.elementIdCounter += self.elementIdCounterIncr
			return nid

	def AddFeature(self, feature, options):
		geometry = feature.geometry
		if isinstance(geometry, geom.Point):
			geometry.id = self.elementIdCounter
			self.elementIdCounter += self.elementIdCounterIncr
			self.features.append(feature)

		if isinstance(geometry, geom.Way):
			geometry.id = self.elementIdCounter
			assert geometry.id != 0
			self.elementIdCounter += self.elementIdCounterIncr

			for pt in geometry.points:
				pt.id = self.DetermineNodeId(pt, options)
			self.features.append(feature)

		if isinstance(geometry, geom.Relation):
			geometry.id = self.elementIdCounter
			self.elementIdCounter += self.elementIdCounterIncr

			members = feature.members
			for mem in members:
				print mem

			self.features.append(feature)

	
