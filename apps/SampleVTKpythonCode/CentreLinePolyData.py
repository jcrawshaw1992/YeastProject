def CentreLinePolyData(self):
        """Compute centrelines based on the profile, reusing our
        memoed copy or reading from the cache file if possible.
        """
        if (
            os.path.exists(self.CentreLineFile)
            and os.path.getmtime(self.CentreLineFile) > os.path.getmtime(self.StlFile)
            and os.path.getmtime(self.CentreLineFile) > os.path.getmtime(self.FileName)
        ):
            # Cached!
            reader = vtk.vtkXMLPolyDataReader()
            reader.SetFileName(self.CentreLineFile)
            reader.Update()
            return reader.GetOutput()

        # Have to compute it

        # Read the STL file
        reader = vtk.vtkSTLReader()
        reader.SetFileName(profile.StlFile)

        # Find the seed points for centreline calculation
        # Use points one iolet radius back along the normal.

        outletPts = []

        def scale(iolet):
            pt = iolet.Centre - iolet.Radius * iolet.Normal
            pt = pt / self.LengthUnit
            return pt.magnitude

        for iolet in self.Iolets:
            if isinstance(iolet._iolet, Inlet):
                inletPt = scale(iolet)
            else:
                outletPts.append(scale(iolet))
                pass
            continue

        srcPts, tgtPts = FindSeeds(reader, inletPt, outletPts)

        # Lazy import since it's so slow!
        from vmtk import vtkvmtk

        centreliner = vtkvmtk.vtkvmtkPolyDataCenterlines()
        centreliner.SetInputConnection(reader.GetOutputPort())

        centreliner.SetSourceSeedIds(srcPts)
        centreliner.SetTargetSeedIds(tgtPts)
        centreliner.SetRadiusArrayName("radius")
        centreliner.SetCostFunction("1/R")

        cleaner = vtk.vtkCleanPolyData()
        cleaner.SetInputConnection(centreliner.GetOutputPort())

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetInputConnection(cleaner.GetOutputPort())
        writer.SetFileName(self.CentreLineFile)
        writer.Write()
        return cleaner.GetOutput()