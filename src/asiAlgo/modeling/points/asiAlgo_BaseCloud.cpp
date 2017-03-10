//-----------------------------------------------------------------------------
// Created on: 28 November 2016
// Created by: Quaoar
//-----------------------------------------------------------------------------
// Web: http://dev.opencascade.org/, http://quaoar.su/blog
//-----------------------------------------------------------------------------

// Own include
#include <asiAlgo_BaseCloud.h>

// Instantiate for allowed types
template class asiAlgo_BaseCloud<double>;
template class asiAlgo_BaseCloud<float>;

//-----------------------------------------------------------------------------

template <typename TCoordType>
asiAlgo_BaseCloud<TCoordType>::asiAlgo_BaseCloud() {}

//-----------------------------------------------------------------------------

template <typename TCoordType>
asiAlgo_BaseCloud<TCoordType>::asiAlgo_BaseCloud(const std::vector<TCoordType>& coords)
: m_coords(coords) {}

//-----------------------------------------------------------------------------

template <typename TCoordType>
void asiAlgo_BaseCloud<TCoordType>::CopyTo(asiAlgo_BaseCloud<TCoordType>& copy) const
{
  copy.m_coords = m_coords;
}

//-----------------------------------------------------------------------------

template <typename TCoordType>
void asiAlgo_BaseCloud<TCoordType>::Reserve(const int nElems)
{
  m_coords.resize(nElems*3);
}

//-----------------------------------------------------------------------------

template <typename TCoordType>
int asiAlgo_BaseCloud<TCoordType>::GetNumberOfElements() const
{
  return (int) (m_coords.size() / 3);
}

//-----------------------------------------------------------------------------

template <typename TCoordType>
void asiAlgo_BaseCloud<TCoordType>::AddElement(const TCoordType x,
                                               const TCoordType y,
                                               const TCoordType z)
{
  m_coords.push_back(x);
  m_coords.push_back(y);
  m_coords.push_back(z);
}

//-----------------------------------------------------------------------------

template <typename TCoordType>
void asiAlgo_BaseCloud<TCoordType>::SetElement(const int        elemIndex,
                                               const TCoordType x,
                                               const TCoordType y,
                                               const TCoordType z)
{
  m_coords[elemIndex*3 + 0] = x;
  m_coords[elemIndex*3 + 1] = y;
  m_coords[elemIndex*3 + 2] = z;
}

//-----------------------------------------------------------------------------

template <typename TCoordType>
void asiAlgo_BaseCloud<TCoordType>::GetElement(const int   elemIndex,
                                               TCoordType& x,
                                               TCoordType& y,
                                               TCoordType& z) const
{
  x = m_coords[3*elemIndex + 0];
  y = m_coords[3*elemIndex + 1];
  z = m_coords[3*elemIndex + 2];
}

//-----------------------------------------------------------------------------

template <typename TCoordType>
const std::vector<TCoordType>& asiAlgo_BaseCloud<TCoordType>::GetCoords() const
{
  return m_coords;
}

//-----------------------------------------------------------------------------

template <typename TCoordType>
std::vector<TCoordType>& asiAlgo_BaseCloud<TCoordType>::ChangeCoords()
{
  return m_coords;
}

//-----------------------------------------------------------------------------

template <typename TCoordType>
void asiAlgo_BaseCloud<TCoordType>::ComputeBoundingBox(TCoordType& xMin, TCoordType& xMax,
                                                       TCoordType& yMin, TCoordType& yMax,
                                                       TCoordType& zMin, TCoordType& zMax) const
{
  TCoordType _xMin = std::numeric_limits<TCoordType>::max();
  TCoordType _yMin = _xMin;
  TCoordType _zMin = _xMin;
  TCoordType _xMax = std::numeric_limits<TCoordType>::min();
  TCoordType _yMax = _xMax;
  TCoordType _zMax = _xMax;

  const int nElems = this->GetNumberOfElements();
  //
  for ( int e = 0; e < nElems; ++e )
  {
    TCoordType x, y, z;
    this->GetElement(e, x, y, z);

    if ( x > _xMax )
      _xMax = x;
    if ( x < _xMin )
      _xMin = x;
    if ( y > _yMax )
      _yMax = y;
    if ( y < _yMin )
      _yMin = y;
    if ( z > _zMax )
      _zMax = z;
    if ( z < _zMin )
      _zMin = z;
  }

  xMin = _xMin;
  xMax = _xMax;
  yMin = _yMin;
  yMax = _yMax;
  zMin = _zMin;
  zMax = _zMax;
}

//-----------------------------------------------------------------------------

//! Reads base cloud recorded in the input file with common XYZ format. That
//! is, the file contains just coordinate triples without any additional
//! structuring information.
//! \param filename [in] file to read.
//! \return true in case of success, false -- otherwise.
template <typename TCoordType>
bool asiAlgo_BaseCloud<TCoordType>::Load(const char* filename)
{
  std::ifstream FILE(filename);
  if ( !FILE.is_open() )
    return false;

  while ( !FILE.eof() )
  {
    char str[256];
    FILE.getline(str, 256);

    std::vector<std::string> tokens;
    std::istringstream iss(str);
    std::copy( std::istream_iterator<std::string>(iss),
               std::istream_iterator<std::string>(),
               std::back_inserter< std::vector<std::string> >(tokens) );

    if ( tokens.empty() || tokens.size() < 3 )
      continue;

    TCoordType x = (TCoordType) ( atof(tokens[0].c_str()) ),
               y = (TCoordType) ( atof(tokens[1].c_str()) ),
               z = (TCoordType) ( atof(tokens[2].c_str()) );

    this->AddElement(x, y, z);
  }

  FILE.close();
  return true;
}

//-----------------------------------------------------------------------------

//! Writes base cloud to file with given filename.
//! \param filename [in] file to write into.
//! \return true in case of success, false -- otherwise.
template <typename TCoordType>
bool asiAlgo_BaseCloud<TCoordType>::SaveAs(const char* filename) const
{
  std::ofstream FILE(filename);
  if ( !FILE.is_open() )
    return false;

  for ( int e = 0; e < this->GetNumberOfElements(); ++e )
  {
    TCoordType x, y, z;
    this->GetElement(e, x, y, z);

    FILE << x << " " << y << " " << z << "\n";
  }

  FILE.close();
  return true;
}