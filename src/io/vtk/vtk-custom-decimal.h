#include <xsd/cxx/xml/string.hxx> // xml::transcode

#include <xsd/cxx/tree/text.hxx>  // text_content

using namespace std;

// Define decimal as long double
//
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

        static type
        create (const xercesc::DOMElement& e, flags f, type* c)
        {
          return create (text_content<char> (e), 0, f, c);
        }

        static type
        create (const xercesc::DOMAttr& a, flags f, type* c)
        {
          return create (xml::transcode<char> (a.getValue ()), 0, f, c);
        }

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

namespace XERCES_CPP_NAMESPACE
{
  void
  operator<< (xercesc::DOMElement& e, const xml_schema::as_decimal& d);

  void
  operator<< (xercesc::DOMAttr& a, const xml_schema::as_decimal& d);
}

namespace xsd
{
  namespace cxx
  {
    namespace tree
    {
      void
      operator<< (xml_schema::list_stream& ls, const xml_schema::as_decimal& d);
    }
  }
}
