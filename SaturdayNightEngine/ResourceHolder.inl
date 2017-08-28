using namespace std;

template< typename Resource, typename Identifier >
void ResourceHolder< Resource, Identifier >::load( Identifier id, const string& filename )
{
  unique_ptr< Resource > resource( new Resource() );
  bool result = resource->loadFromFile( filename );
  if( !result )
  {
    throw runtime_error( "Failed to load resource " + filename );
  }

  auto insertion = mResourceMap_.insert( make_pair( id, move( resource ) ) );
  assert( insertion.second );
}

template< typename Resource, typename Identifier >
template< typename Parameter >
void ResourceHolder< Resource, Identifier >::load( Identifier id, const string& filename, const Parameter& parameter )
{
  unique_ptr< Resource > resource( new Resource() );
  bool result = texture->loadFromFile( filename, parameter );
  if( !result )
  {
    throw runtime_error( "Failed to load resource " + filename );
  }

  auto insertion = mResourceMap_.insert( make_pair( id, move( texture ) ) );
  assert( insertion.second );
}

template< typename Resource, typename Identifier >
Resource& ResourceHolder< Resource, Identifier >::get( Identifier id )
{ 
  auto texture = mResourceMap_.find( id );
  assert( texture != mResourceMap_.end() );
  return *texture->second;
}

template< typename Resource, typename Identifier >
const Resource& ResourceHolder< Resource, Identifier >::get( Identifier id ) const
{ 
  auto texture = mResourceMap_.find( id );
  assert( texture != mResourceMap_.end() );
  return *texture->second;
}