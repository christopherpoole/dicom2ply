# License & Copyright
# ===================
#
# Copyright 2012 Christopher M Poole
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
# Author:    Christopher M Poole
# Email:     mail@christopherpoole.net


import os
import sys

import dicom
import numpy
import pylab

from PIL import Image, ImageDraw
   

class Contour(object):
    def __init__(self, contour, dicom_dir='', bins=2**12):
        self.dicom_dir = dicom_dir
        self.bins = bins
        self.slice_ref = contour.ContourImages[0].RefdSOPInstanceUID
                
        self.y = contour.ContourData[0::3]
        self.x = contour.ContourData[1::3]
        self.z = contour.ContourData[2::3]
        
        self.vertex_count = len(self.x)
        
        self._stats()        
        
    def _stats(self):
        self.slice_file = dicom.read_file("%s/CT.%s.dcm" % (self.dicom_dir, self.slice_ref))
        self.slice_data = self.slice_file.PixelArray
        
        x = numpy.array(self.x)# + 256
        y = numpy.array(self.y) + self.slice_file.TableHeight
                
        self.mask = self._get_mask(x, y)
        
        self.masked = self.slice_data * self.mask
        self.masked = self.masked.flatten()
        self.masked = self.masked[self.masked.nonzero()[0]]
        
        if self.masked.size == 0:
            self.histogram = None
            self.mode = None
            self.mean = None
            self.std = None
            self.median = None
        else:
            self.histogram = numpy.histogram(self.masked, bins=self.bins)
            self.mode = self.histogram[1][numpy.argmax(self.histogram[0])]
            self.mean = numpy.mean(self.masked)
            self.std = numpy.std(self.masked)
            self.median = numpy.median(self.masked) 
        
    def _get_mask(self, x, y, shape=(512, 512)):
        pylab.plot(x, y)
    
        mask = numpy.zeros(shape)
        points = []
        for i in range(0, len(x)):
            points.append((x[i], y[i]))
        
        if len(points) < 3:
            return mask
        
        im = Image.fromarray(mask)
        draw = ImageDraw.Draw(im)
        draw.polygon(points, fill=1)
        mask = numpy.asarray(im)
        
        return mask.astype(numpy.int8)
        
        
class RegionOfInterest(object):
    def __init__(self, roi, name, dicom_dir='', bins=2**12):
        self.dicom_dir = dicom_dir
        self.bins = bins
                
        self.name = name
                
        self.contours = []
        self.vertex_count = 0
        for contour in roi.Contours:
            c = Contour(contour, dicom_dir=self.dicom_dir)
            if c.mean is None:
                continue
            self.contours.append(c)
            self.vertex_count += c.vertex_count       
        
        self._stats()
    
    
    def __len__(self):
        return len(self.contours)
    
    
    def _stats(self):
        masked = numpy.empty((0))
        
        for c in self.contours:
            masked = numpy.append(masked, c.masked)
        
        try:
            self.histogram = numpy.histogram(masked.flatten(), bins=self.bins)
            self.mode = self.histogram[1][numpy.argmax(self.histogram[0])]
            self.mean = numpy.mean(masked)
            self.std = numpy.std(masked)
            self.median = numpy.median(masked)
            self.sum = numpy.sum(masked)
            self.len = len(masked)
        except ValueError:
            return
            
    @property
    def mask_sum(self):
        total = 0
        for c in self.contours:
            total += c.mask.sum()
        return total
        
        
    @property
    def masked_sum(self):
        total = 0
        for c in self.contours:
            total += c.masked.sum()
        return total 


class Patient(object):
    def __init__(self, dicom_dir, debug=True):
        self.debug = debug
        self.dicom_dir = dicom_dir
        
        self.files = os.walk(self.dicom_dir).next()[2]
            
        for f in self.files:
            if f[:2] == 'RS':
                dicom_structure = f
                print f
                break
                
        self.structure = dicom.ReadFile("%s/%s" % (dicom_dir, dicom_structure))
        
        self.region_names = {}
        for i, roi in enumerate(self.structure.RTROIObservations):
            try:
                name = roi.ROIObservationLabel
            except AttributeError:
                name = i
            number = roi.ObservationNumber
            self.region_names[number] = name
        
        self.regions = {}
        for roi in self.structure.ROIContours:
            number = roi.ReferencedROINumber
            try:
                name = self.region_names[number]
            except KeyError:
                # ROI did not have a label
                continue
            
            try:
                roi.Contours
            except AttributeError:
                # Not an ROI with geometry
                continue

            r = RegionOfInterest(roi, name, dicom_dir)
            self.regions[self.region_names[number]] = r

   
    def dump_ply(self, directory='.', names=[]):
        if names == []:
            names = self.regions.keys()
          
        for n in names:
            print n   
            roi = self.regions[n]
            
            ply_verts = []
            for contour in roi.contours:
                for i in range(contour.vertex_count):
                    ply_verts.append("%f %f %f" % (contour.x[i], contour.y[i], contour.z[i]))
            
            try:
                ply_header = ['ply',
                    'format ascii 1.0',
                    'comment name roi_%s' % roi.name,
                    'comment mean %f' % roi.mean,
                    'comment std %f' % roi.std,
                    'comment median %f' % roi.median,
                    'comment mode %f' % roi.mode,
                    'comment sum %f' % roi.sum,
                    'comment len %f' % roi.len,
                    'element vertex %i' % roi.vertex_count,
                    'property float x',
                    'property float y',
                    'property float z',
                    'end_header']
            except AttributeError:
                print " - fail"
                continue           
            file_name = "%s/roi_%s.ply" % (directory, roi.name)
            print file_name
            ply_file = file(file_name, 'w')
            for line in ply_header:
                ply_file.write("%s\n" % line)
            for line in ply_verts:
                ply_file.write("%s\n" % line)    
            ply_file.close()
    
    @property
    def roi_names(self):
        return self._rois.keys() 
    
    
if __name__ == "__main__":
    patient = Patient(sys.argv[1])
    patient.dump_ply(directory=sys.argv[2])

