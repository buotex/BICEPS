//
// XMLWriter.cpp
//
//
// Original author: Darren Kessner <Darren.Kessner@cshs.org>
//
// Copyright 2007 Spielberg Family Center for Applied Proteomics
//   Cedars-Sinai Medical Center, Los Angeles, California  90048
//
// Licensed under the Apache License, Version 2.0 (the "License"); 
// you may not use this file except in compliance with the License. 
// You may obtain a copy of the License at 
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software 
// distributed under the License is distributed on an "AS IS" BASIS, 
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
// See the License for the specific language governing permissions and 
// limitations under the License.
//

#define PWIZ_SOURCE

#include "XMLWriter.hpp"
#include "utility/misc/Stream.hpp"
#include "utility/misc/String.hpp"
#include "utility/misc/Exception.hpp"
#include <stack>
using std::stack;

namespace pwiz {
namespace minimxml {


class XMLWriter::Impl
{
    public:

    Impl(ostream& os, const Config& config);
    void pushStyle(unsigned int flags);
    void popStyle(); 
    void processingInstruction(const string& name, const string& data);
    void startElement(const string& name, 
                      const Attributes& attributes,
                      EmptyElementTag emptyElementTag);
    void endElement();
    void characters(const string& text);
    bio::stream_offset position() const;
    bio::stream_offset positionNext() const;

    private:
    ostream& os_;
    Config config_;
    stack<string> elementStack_;
    stack<unsigned int> styleStack_;

    string indentation() const {return string(elementStack_.size()*config_.indentationStep, ' ');}
    string indentation(size_t depth) const {return string(depth*config_.indentationStep, ' ');}
    bool style(StyleFlag styleFlag) const {return styleStack_.top() & styleFlag ? true : false;}
};


XMLWriter::Impl::Impl(ostream& os, const Config& config)
:   os_(os), config_(config)
{
    styleStack_.push(config.initialStyle);
}


void XMLWriter::Impl::pushStyle(unsigned int flags)
{
    styleStack_.push(flags);
}


void XMLWriter::Impl::popStyle() 
{
    styleStack_.pop();
    if (styleStack_.empty())
        throw runtime_error("[XMLWriter] Style stack underflow.");
}


void XMLWriter::Impl::processingInstruction(const string& name, const string& data) 
{
    ostream* os = &os_;
    if (config_.outputObserver) os = new ostringstream;

    *os << indentation() << "<?" << name << " " << data << "?>\n"; 

    if (config_.outputObserver)
    {
        config_.outputObserver->update(static_cast<ostringstream*>(os)->str());
        os_ << static_cast<ostringstream*>(os)->str();
        delete os;
    }
}


void writeEscapedAttributeXML(ostream& os, const string& str)
{
    for (size_t i=0; i < str.size(); ++i)
    {
        const char& c = str[i];
        switch (c)
        {
            case '&': os << "&amp;"; break;
            case '"': os << "&quot;"; break;
            case '\'': os << "&apos;"; break;
            case '<': os << "&lt;"; break;
            case '>': os << "&gt;"; break;
            default: os << c; break;
        }
    }
}


void writeEscapedTextXML(ostream& os, const string& str)
{
    for (size_t i=0; i < str.size(); ++i)
    {
        const char& c = str[i];
        switch (c)
        {
            case '&': os << "&amp;"; break;
            case '<': os << "&lt;"; break;
            case '>': os << "&gt;"; break;
            default: os << c; break;
        }
    }
}


void XMLWriter::Impl::startElement(const string& name, 
                  const Attributes& attributes,
                  EmptyElementTag emptyElementTag)
{
    ostream* os = &os_;
    if (config_.outputObserver) os = new ostringstream;

    if (!style(StyleFlag_InlineOuter))
        *os << indentation();

    *os << "<" << name;

    string attributeIndentation(name.size()+1, ' ');

    for (Attributes::const_iterator it=attributes.begin(); it!=attributes.end(); ++it)
    {
        *os << " " << it->first << "=\"";
        writeEscapedAttributeXML(*os, it->second);
        *os << "\"";
        if (style(StyleFlag_AttributesOnMultipleLines) && (it+1)!=attributes.end())
            *os << "\n" << indentation() << attributeIndentation;
    }

    *os << (emptyElementTag==EmptyElement ? "/>" : ">");

    if (!style(StyleFlag_InlineInner) || 
        !style(StyleFlag_InlineOuter) && emptyElementTag==EmptyElement)
        *os << "\n";

    if (emptyElementTag == NotEmptyElement)
        elementStack_.push(name);

    if (config_.outputObserver)
    {
        config_.outputObserver->update(static_cast<ostringstream*>(os)->str());
        os_ << static_cast<ostringstream*>(os)->str();
        delete os;
    }
}


void XMLWriter::Impl::endElement()
{
    ostream* os = &os_;
    if (config_.outputObserver) os = new ostringstream;

    if (elementStack_.empty())
        throw runtime_error("[XMLWriter] Element stack underflow.");

    if (!style(StyleFlag_InlineInner))
        *os << indentation(elementStack_.size()-1);

    *os << "</" << elementStack_.top() << ">";
    elementStack_.pop();

    if (!style(StyleFlag_InlineOuter))
        *os << "\n";
        
    if (config_.outputObserver)
    {
        config_.outputObserver->update(static_cast<ostringstream*>(os)->str());
        os_ << static_cast<ostringstream*>(os)->str();
        delete os;
    }
}


void XMLWriter::Impl::characters(const string& text)
{
    ostream* os = &os_;
    if (config_.outputObserver) os = new ostringstream;

    if (!style(StyleFlag_InlineInner))
        *os << indentation();

    writeEscapedTextXML(*os, text);

    if (!style(StyleFlag_InlineInner))
        *os << "\n";

    if (config_.outputObserver)
    {
        config_.outputObserver->update(static_cast<ostringstream*>(os)->str());
        os_ << static_cast<ostringstream*>(os)->str();
        delete os;
    }
}


XMLWriter::stream_offset XMLWriter::Impl::position() const
{
    os_ << flush;
    return boost::iostreams::position_to_offset(os_.tellp()); 
}


XMLWriter::stream_offset XMLWriter::Impl::positionNext() const
{
    stream_offset offset = position(); 
    if (!style(StyleFlag_InlineOuter))
        offset += indentation().size();
    return offset;
}


//
// XMLWriter forwarding functions 
//


PWIZ_API_DECL XMLWriter::XMLWriter(ostream& os, const Config& config)
:   impl_(new Impl(os, config))
{}

PWIZ_API_DECL void XMLWriter::pushStyle(unsigned int flags) {impl_->pushStyle(flags);}

PWIZ_API_DECL void XMLWriter::popStyle() {impl_->popStyle();}

PWIZ_API_DECL void XMLWriter::processingInstruction(const string& name, const string& data) 
{
    impl_->processingInstruction(name, data);
}

PWIZ_API_DECL void XMLWriter::startElement(const string& name, 
                             const Attributes& attributes,
                             EmptyElementTag emptyElementTag)
{
    impl_->startElement(name, attributes, emptyElementTag);
}

PWIZ_API_DECL void XMLWriter::endElement() {impl_->endElement();}

PWIZ_API_DECL void XMLWriter::characters(const string& text) {impl_->characters(text);}

PWIZ_API_DECL XMLWriter::stream_offset XMLWriter::position() const {return impl_->position();}

PWIZ_API_DECL XMLWriter::stream_offset XMLWriter::positionNext() const {return impl_->positionNext();}


} // namespace minimxml
} // namespace pwiz


