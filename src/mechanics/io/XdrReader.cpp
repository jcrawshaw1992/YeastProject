// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "xdr.hpp"
#include "XdrReader.hpp"
#include "Debug.hpp"

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      namespace xdr
      {

        XdrReader::XdrReader()
        {
          TRACE("XdrReader constructor here")
        }

        // Functions to read out the next bit of the file as a certain type.
        bool XdrReader::readDouble(double& outDouble)
        {
          TRACE("readDouble")
          return xdr_double(&mXdr, &outDouble);
        }

        bool XdrReader::readFloat(float& outDouble)
        {
          TRACE("readFloat")
          return xdr_float(&mXdr, &outDouble);
        }

        bool XdrReader::readInt(int& outInt)
        {
           TRACE("readInt")
          return xdr_int(&mXdr, &outInt);
        }

        bool XdrReader::readUnsignedInt(unsigned int& outUInt)
        {
          TRACE("readUnsignedInt")
          return xdr_u_int(&mXdr, &outUInt);
        }

        bool XdrReader::readString(std::string& outString, unsigned maxLength)
        {
          TRACE("readString")
          char* temp = new char[maxLength];
          bool ret = xdr_string(&mXdr, &temp, maxLength);
          outString = std::string(temp);
          delete temp;
          return ret;
        }

        bool XdrReader::readUnsignedLong(uint64_t& outULong)
        {
          TRACE("Jess got here. good jess")
          u_quad_t temporary;
#ifdef __APPLE__
          bool ret = xdr_u_int64_t(&mXdr, &temporary);
#else
          bool ret = xdr_uint64_t(&mXdr, &temporary);
#endif
          outULong = temporary;
          return ret;
        }

        unsigned int XdrReader::GetPosition()
        {
          return xdr_getpos(&mXdr);
        }

        // Returns false on failure
        bool XdrReader::SetPosition(unsigned int iPosition)
        {
          return xdr_setpos(&mXdr, iPosition);
        }

        // Destructor to get rid of any resources used by the Xdr object. This class doesn't create
        // the file object, so it doesn't free it either.
        XdrReader::~XdrReader()
        {
          xdr_destroy(&mXdr);
        }

      } // namespace xdr
    } // namespace writers
  }
}
