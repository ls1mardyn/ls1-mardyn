/**
 * @file vtk-custom-decimal.h
 * @brief Required for mapping the 'decimal' data type to 'long double' (would be 'double' otherwise).
 */

#include <xsd/cxx/xml/string.hxx> // xml::transcode

#include <xsd/cxx/tree/text.hxx>  // text_content

// Define decimal as long double
//
/**
 * @brief Include 'long double' as valid data types 
 * 
 */
namespace xsd
{
  namespace cxx
  {
    namespace tree
    {
      template <>
      struct fundamental_p<long double>
      {
        static const bool r = true;
      };
    }
  }
}

// Parsing
//
/**
 * @brief Enables the parsing of XML data to the 'long double' type.
 * 
 */
namespace xsd
{
  namespace cxx
  {
    namespace tree
    {
      
      template <>
      struct traits<long double, char, schema_type::decimal>
      {
        
        typedef long double type;

        /**
         * @brief Creates a long double instance from a DOM element.
         * 
         * @param e A DOM element to extract the data from.
         * @param f Flags to create the new instance with.
         * @param c A pointer to the object that will contain the new instance.
         * @return type = long double
         */
        static type
        create (const xercesc::DOMElement& e, flags f, type* c)
        {
          return create (text_content<char> (e), 0, f, c);
        }

        /**
         * @brief Creates a long double instance from a DOM attribute.
         * 
         * @param a A DOM attribute to extract the data from.
         * @param f Flags to create the new instance with.
         * @param c A pointer to the object that will contain the new instance.
         * @return type = long double
         */
        static type
        create (const xercesc::DOMAttr& a, flags f, type* c)
        {
          return create (xml::transcode<char> (a.getValue ()), 0, f, c);
        }
        
        /**
         * @brief Creates a long double instance from a string fragment.
         * 
         * @param s A string fragment to extract the data from.
         * @return type = long double 
         */
        static type
        create (const std::string& s,
                const xercesc::DOMElement*,
                flags,
                type*);
      };
    }
  }
}

// Serialization
//
/**
 * @brief Enables the serialization of the 'long double' type to XML.
 * 
 */

namespace XERCES_CPP_NAMESPACE
{
  /**
   * @brief Insertion operator of a long double to a DOM element.
   * 
   * @param e Reference to a DOM element to serialize the data to.
   * @param d Reference to a decimal (= long double) supposed to be serialized.
   */
  void
  operator<< (xercesc::DOMElement& e, const xml_schema::as_decimal& d);

  /**
   * @brief Insertion operator of a long double to a DOM attribute.
   * 
   * @param a Reference to a DOM attribute to serialize the data to.
   * @param d Reference to a decimal (= long double) supposed to be serialized.
   */
  void
  operator<< (xercesc::DOMAttr& a, const xml_schema::as_decimal& d);
}

namespace xsd
{
  namespace cxx
  {
    namespace tree
    {
      /**
       * @brief Insertion operator of a long double to a list stream (wrapper class around basic_ostringstream).
       * 
       * @param ls Reference to list_stream to serialize the data to.
       * @param d Reference to a decimal (= long double) supposed to be serialized.
       */
      void
      operator<< (xml_schema::list_stream& ls, const xml_schema::as_decimal& d);
    }
  }
}
