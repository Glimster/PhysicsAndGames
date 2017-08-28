#pragma once

template< typename Resource, typename Identifier >
class ResourceHolder
{
public:
  void load( Identifier id, const std::string& filename );
  
  // Useful when passing parameter to LoadFromFile, e.g. SFML shader
  template< typename Parameter >
  void load( Identifier id, const std::string& filename, const Parameter& parameter );

  Resource& get( Identifier id );
  const Resource& get( Identifier id ) const;

private:
  std::map< Identifier, std::unique_ptr< Resource > > mResourceMap_;
};

#include "ResourceHolder.inl"