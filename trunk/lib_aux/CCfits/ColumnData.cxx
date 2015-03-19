//	Astrophysics Science Division,
//	NASA/ Goddard Space Flight Center
//	HEASARC
//	http://heasarc.gsfc.nasa.gov
//	e-mail: ccfits@legacy.gsfc.nasa.gov
//
//	Original author: Ben Dorman
#ifdef _MSC_VER
#include "MSconfig.h" // for truncation warning
#endif

// ColumnData
#include "ColumnData.h"
#include <algorithm>
namespace CCfits
{

#ifndef SPEC_TEMPLATE_DECL_DEFECT
        template <>
        void ColumnData<String>::writeData (const std::vector<String>& indata, 
                        long firstRow, String* nullValue)
        {
          int    status=0;
          char** columnData=FITSUtil::CharArray(indata);

          if ( fits_write_colnull(fitsPointer(), TSTRING, index(), firstRow, 1, indata.size(),
			          columnData, 0, &status) != 0 )
              throw FitsError(status);
          unsigned long elementsToWrite (indata.size() + firstRow - 1);
          std::vector<String> __tmp(m_data);
          if (m_data.size() < elementsToWrite) 
          {
                  m_data.resize(elementsToWrite,"");
                  std::copy(__tmp.begin(),__tmp.end(),m_data.begin());
          }
          std::copy(indata.begin(),indata.end(),m_data.begin()+firstRow-1);


          for (size_t i = 0; i < indata.size(); ++i)
          {
                delete [] columnData[i];
          }
          delete [] columnData;
        }  

#endif

#ifndef SPEC_TEMPLATE_IMP_DEFECT
#ifndef SPEC_TEMPLATE_DECL_DEFECT
        template <>
        void ColumnData<String>::readColumnData (long firstRow, 
                                        long nelements, 
                                        String* nullValue)
        {
          int status = 0;

           int   anynul = 0;
           char** array = new char*[nelements]; 

           int j(0);
           for ( ; j < nelements; ++j)
           {
               array[j] = new char[width() + 1];
           }

           char* nulval = 0;
           if (nullValue) 
           {
                   nulval = const_cast<char*>(nullValue->c_str());
           }
           else
           {
                nulval = new char;
                *nulval = '\0';       
           }


          try
          {
                makeHDUCurrent();
                if (fits_read_col_str(fitsPointer(),index(), firstRow,1,nelements,
                  nulval,array, &anynul,&status) ) throw FitsError(status);
          }
          catch (FitsError)
          {
                // ugly. but better than leaking resources.       
                for (int jj = 0; jj < nelements; ++jj)
                {
                        delete [] array[jj];
                }     

                delete [] array; 
                delete nulval;
                throw; 
          }


          if (m_data.size() != rows()) setData(std::vector<String>(rows(),String(nulval)));

          // the 'first -1 ' converts to zero based indexing.

          for ( j = 0; j < nelements; j++)
          {
                m_data[j - 1 + firstRow] = String(array[j]);
          }

          for ( j = 0; j < nelements; j++)
          {
                delete [] array[j];
          }     

          delete [] array; 
          delete nulval; 
          if (nelements == rows()) isRead(true); 

        }
#endif
#endif

#ifndef SPEC_TEMPLATE_IMP_DEFECT
#ifndef SPEC_TEMPLATE_DECL_DEFECT
template <>
void ColumnData<complex<float> >::setDataLimits (complex<float>* limits)
       {
                m_minLegalValue = limits[0];
                m_maxLegalValue = limits[1];
                m_minDataValue =  limits[2];
                m_maxDataValue =  limits[3];
        }
template <>
void ColumnData<complex<double> >::setDataLimits (complex<double>* limits)
        {
                m_minLegalValue = limits[0];
                m_maxLegalValue = limits[1];
                m_minDataValue =  limits[2];
                m_maxDataValue =  limits[3];
        }


        template <>
        void ColumnData<complex<float> >::readColumnData (long firstRow,
                                                long nelements,
                                                complex<float>* nullValue)
        {
          // specialization for ColumnData<String> 
          int status(0);
          int   anynul(0);
          FITSUtil::auto_array_ptr<float> pArray(new float[nelements*2]); 
          float* array = pArray.get();
          float nulval(0);
          makeHDUCurrent();


          if (fits_read_col_cmp(fitsPointer(),index(), firstRow,1,nelements,
                  nulval,array, &anynul,&status) ) throw FitsError(status);


          if (m_data.size() != rows()) m_data.resize(rows());

          // the 'j -1 ' converts to zero based indexing.

          for (int j = 0; j < nelements; ++j)
          {

                m_data[j - 1 +firstRow] = std::complex<float>(array[2*j],array[2*j+1]);
          }

          if (nelements == rows()) isRead(true); 
        }

        template <>
        void ColumnData<complex<double> >::readColumnData (long firstRow, 
                                                        long nelements,
                                                        complex<double>* nullValue)
        {
          // specialization for ColumnData<complex<double> > 
           int status(0);
           int   anynul(0);
           FITSUtil::auto_array_ptr<double> pArray(new double[nelements*2]); 
           double* array = pArray.get();
           double nulval(0);
           makeHDUCurrent();


          if (fits_read_col_dblcmp(fitsPointer(), index(), firstRow,1,nelements,
                  nulval,array, &anynul,&status) ) throw FitsError(status);




          if (m_data.size() != rows()) setData(std::vector<complex<double> >(rows(),nulval));

          // the 'j -1 ' converts to zero based indexing.

          for (int j = 0; j < nelements; j++)
          {

                m_data[j - 1 + firstRow] = std::complex<double>(array[2*j],array[2*j+1]);
          }
          if (nelements == rows()) isRead(true); 

        }
#endif
#endif

#ifndef SPEC_TEMPLATE_DECL_DEFECT
        template <>
        void ColumnData<complex<float> >::writeData (const std::vector<complex<float> >& inData, 
                                long firstRow, 
                                complex<float>* nullValue)
        {
                int status(0);
                int nRows (inData.size());
                FITSUtil::auto_array_ptr<float> pData(new float[nRows*2]);
                float* Data = pData.get();
                std::vector<complex<float> > __tmp(m_data);
                for (int j = 0; j < nRows; ++j)
                {
                        Data[ 2*j] = inData[j].real();
                        Data[ 2*j + 1] = inData[j].imag();
                }     

                try
                {

                        if (fits_write_col_cmp(fitsPointer(), index(), firstRow, 1, 
                                  nRows,Data, &status) != 0) throw FitsError(status);
                        long elementsToWrite(nRows + firstRow -1);
                        if (elementsToWrite > static_cast<long>(m_data.size())) 
                        {

                                m_data.resize(elementsToWrite);
                        }

                        std::copy(inData.begin(),inData.end(),m_data.begin()+firstRow-1);

                        // tell the Table that the number of rows has changed
                        parent()->updateRows();
                }
                catch (FitsError) // the only thing that can throw here.
                {
                        // reset to original content and rethrow the exception.
                        m_data.resize(__tmp.size());
                        m_data = __tmp;
                }      

        }
#endif

#ifndef SPEC_TEMPLATE_DECL_DEFECT
        template <>
        void ColumnData<complex<double> >::writeData (const std::vector<complex<double> >& inData, 
                                long firstRow, 
                                complex<double>* nullValue)
        {
                int status(0);
                int nRows (inData.size());
                FITSUtil::auto_array_ptr<double> pData(new double[nRows*2]);
                double* Data = pData.get();
                std::vector<complex<double> > __tmp(m_data);
                for (int j = 0; j < nRows; ++j)
                {
                        pData[ 2*j] = inData[j].real();
                        pData[ 2*j + 1] = inData[j].imag();
                }     

                try
                {

                        if (fits_write_col_dblcmp(fitsPointer(), index(), firstRow, 1, 
                                  nRows,Data, &status) != 0) throw FitsError(status);
                        long elementsToWrite(nRows + firstRow -1);
                        if (elementsToWrite > static_cast<long>(m_data.size())) 
                        {

                                m_data.resize(elementsToWrite);
                        }

                        std::copy(inData.begin(),inData.end(),m_data.begin()+firstRow-1);

                        // tell the Table that the number of rows has changed
                        parent()->updateRows();
                }
                catch (FitsError) // the only thing that can throw here.
                {
                        // reset to original content and rethrow the exception.
                        m_data.resize(__tmp.size());
                        m_data = __tmp;
                }      

        }
#endif       

} // namespace CCfits

