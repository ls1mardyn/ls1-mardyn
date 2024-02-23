/** \file Expression.cpp
  * \brief Expression (tree) with variables
  * \author Martin Bernreuther <bernreuther@hlrs.de>
*/

#define EXPRESSION_CPP

#include "Expression.h"



bool Expression::Value::operator==(Value const& v) const
{
	if(isFloat() && v.isFloat())
		return getValueFloat()==v.getValueFloat();
	else if(isFloat() && v.isInt())
		return getValueFloat()==v.getValueInt();
	else if(isInt() && v.isFloat())
		return getValueInt()==v.getValueFloat();
	else if(isInt() && v.isInt())
		return getValueInt()==v.getValueInt();
	return false;
}

bool Expression::Value::operator>(Value const& v) const
{
	if(isFloat() && v.isFloat())
		return getValueFloat()>v.getValueFloat();
	else if(isFloat() && v.isInt())
		return getValueFloat()>v.getValueInt();
	else if(isInt() && v.isFloat())
		return getValueInt()>v.getValueFloat();
	else if(isInt() && v.isInt())
		return getValueInt()>v.getValueInt();
	return false;
}

bool Expression::Value::operator>=(Value const& v) const
{
	if(isFloat() && v.isFloat())
		return getValueFloat()>=v.getValueFloat();
	else if(isFloat() && v.isInt())
		return getValueFloat()>=v.getValueInt();
	else if(isInt() && v.isFloat())
		return getValueInt()>=v.getValueFloat();
	else if(isInt() && v.isInt())
		return getValueInt()>=v.getValueInt();
	return false;
}

bool Expression::Value::operator<(Value const& v) const
{
	if(isFloat() && v.isFloat())
		return getValueFloat()<v.getValueFloat();
	else if(isFloat() && v.isInt())
		return getValueFloat()<v.getValueInt();
	else if(isInt() && v.isFloat())
		return getValueInt()<v.getValueFloat();
	else if(isInt() && v.isInt())
		return getValueInt()<v.getValueInt();
	return false;
}

bool Expression::Value::operator<=(Value const& v) const
{
	if(isFloat() && v.isFloat())
		return getValueFloat()<=v.getValueFloat();
	else if(isFloat() && v.isInt())
		return getValueFloat()<=v.getValueInt();
	else if(isInt() && v.isFloat())
		return getValueInt()<=v.getValueFloat();
	else if(isInt() && v.isInt())
		return getValueInt()<=v.getValueInt();
	return false;
}

const Expression::Value Expression::Value::operator+(Value const& v) const
{
	if(isFloat() && v.isFloat())
		return Value(getValueFloat()+v.getValueFloat());
	else if(isFloat() && v.isInt())
		return Value(getValueFloat()+v.getValueInt());
	else if(isInt() && v.isFloat())
		return Value(getValueInt()+v.getValueFloat());
	else if(isInt() && v.isInt())
		return Value(getValueInt()+v.getValueInt());
	return Value();
}

const Expression::Value Expression::Value::operator-(Value const& v) const
{
	if(isFloat() && v.isFloat())
		return Value(getValueFloat()-v.getValueFloat());
	else if(isFloat() && v.isInt())
		return Value(getValueFloat()-v.getValueInt());
	else if(isInt() && v.isFloat())
		return Value(getValueInt()-v.getValueFloat());
	else if(isInt() && v.isInt())
		return Value(getValueInt()-v.getValueInt());
	return Value();
}

const Expression::Value Expression::Value::operator*(Value const& v) const
{
	if(isFloat() && v.isFloat())
		return Value(getValueFloat()*v.getValueFloat());
	else if(isFloat() && v.isInt())
		return Value(getValueFloat()*v.getValueInt());
	else if(isInt() && v.isFloat())
		return Value(getValueInt()*v.getValueFloat());
	else if(isInt() && v.isInt())
		return Value(getValueInt()*v.getValueInt());
	return Value();
}

const Expression::Value Expression::Value::operator/(Value const& v) const
{
	if(isFloat() && v.isFloat())
		return Value(getValueFloat()/v.getValueFloat());
	else if(isFloat() && v.isInt())
		return Value(getValueFloat()/v.getValueInt());
	else if(isInt() && v.isFloat())
		return Value(getValueInt()/v.getValueFloat());
	else if(isInt() && v.isInt())
		return Value(getValueInt()/v.getValueInt());
	return Value();
}


Expression::Variable* Expression::VariableSet::addVariable(const std::string& name)
{
	std::string vgrpname,varname(name);
	size_t colonpos=name.find(":");
	if(colonpos!=std::string::npos)
	{
		vgrpname=name.substr(0,colonpos);
		varname=name.substr(colonpos+1);
	}
	if(!existVariableGroup(vgrpname))
	{ // create variable group
		//_vargroups[vgrpname]=VariableGroup(vgrpname);
		_vargroups.insert(std::pair<std::string,VariableGroup>(vgrpname,VariableGroup(vgrpname)));
	}
	//_variables[name]=Variable(varname,&_vargroups[vgrpname],val);
	_variables.insert(std::pair<std::string,Variable>(name,Variable(varname,&_vargroups[vgrpname])));
	return &_variables[name];
}

bool Expression::VariableSet::removeVariable(const std::string& name)
{
	Variable* var=getVariable(name);
	if(var)
	{
		const VariableGroup* vargrp=var->getVariableGroup();
		_variables.erase(var->getfullName());
		if(vargrp)
		{
			//vargrp->removeVariable(name);	// done by Variable destructor
			/* erase empty VariableGroup */
			if(!vargrp->countVariables()) _vargroups.erase(vargrp->getName());
		}
		return true;
	}
	else
		return false;
}



void Expression::Node::traverse(std::list<const Node*>& nodelist, enum Etraversetype traversetype) const
{
	switch(traversetype)
	{
		case traversetypePREFIX:
			nodelist.push_back(this);
			if(_children[0]) _children[0]->traverse(nodelist,traversetype);
			if(_children[1]) _children[1]->traverse(nodelist,traversetype);
			break;
		case traversetypeINFIX:
			if(_children[0]) _children[0]->traverse(nodelist,traversetype);
			nodelist.push_back(this);
			if(_children[1]) _children[1]->traverse(nodelist,traversetype);
			break;
		case traversetypePOSTFIX:
			if(_children[0]) _children[0]->traverse(nodelist,traversetype);
			if(_children[1]) _children[1]->traverse(nodelist,traversetype);
			nodelist.push_back(this);
			break;
	}
}

void Expression::Node::writeSubExpr(std::ostream& ostrm, enum Etraversetype traversetype, char sep) const
{
	switch(traversetype)
	{
		case traversetypePREFIX:
			write(ostrm);
			if(sep) ostrm << sep;
			if(_children[0]) _children[0]->writeSubExpr(ostrm,traversetype,sep);
			if(_children[1]) _children[1]->writeSubExpr(ostrm,traversetype,sep);
			break;
		case traversetypeINFIX:
			if(dynamic_cast<const NodeFunction*>(this) != NULL)
			{	// functions format is FUNCTIONNAME(ARGUMENT1, [ARGUMENT2])
				write(ostrm);
				ostrm<<"(";
				if(_children[1]) _children[1]->writeSubExpr(ostrm,traversetype,sep);
				if(_children[0])
				{
					ostrm<<",";
					_children[0]->writeSubExpr(ostrm,traversetype,sep);
				}
				ostrm<<")";
			} else {
				if(_children[0])
				{
					if(_children[0]->_priority>_priority) ostrm<<"(";
					_children[0]->writeSubExpr(ostrm,traversetype,sep);
					if(_children[0]->_priority>_priority) ostrm<<")";
				}
				write(ostrm);
				if(sep) ostrm << sep;
				if(_children[1])
				{
					if(_children[1]->_priority>_priority) ostrm<<"(";
					_children[1]->writeSubExpr(ostrm,traversetype,sep);
					if(_children[1]->_priority>_priority) ostrm<<")";
				}
			}
			break;
		case traversetypePOSTFIX:
			if(_children[0]) _children[0]->writeSubExpr(ostrm,traversetype,sep);
			if(_children[1]) _children[1]->writeSubExpr(ostrm,traversetype,sep);
			write(ostrm);
			if(sep) ostrm << sep;
			break;
	}
}


Expression::NodeOperation2::NodeOperation2(char op, Node* child1, Node* child2, Node* parent)
	: Node(child1,child2,parent), _operator(op)
{
	switch(op)
	{
		case '*':
		case '/':
			_priority=2;
			break;
		case '-':
			_priority=3;
			break;
		case '+':
			_priority=4;
			break;
	}
}

Expression::Value Expression::NodeOperation2::evaluate() const
{
	if(_children[0] && _children[1])
	{
		switch(_operator)
		{
			case '+': return _children[0]->evaluate()+_children[1]->evaluate();
			case '-': return _children[0]->evaluate()-_children[1]->evaluate();
			case '*': return _children[0]->evaluate()*_children[1]->evaluate();
			case '/': return _children[0]->evaluate()/_children[1]->evaluate();
			// the modulo operator % will be realized with a function MOD
			// the power operator ^ (**) will be realized with a function POW
		}
	}
	return Value();
}


Expression::NodeFunction::Efunctype Expression::NodeFunction::functype(const std::string& name)
{
	if(name.compare("abs")==0||name.compare("ABS")==0)
		return functypeABS;
	else if(name.compare("float")==0||name.compare("FLOAT")==0)
		return functypeFLOAT;
	else if(name.compare("int")==0||name.compare("INT")==0)
		return functypeINT;
	else if(name.compare("floor")==0||name.compare("FLOOR")==0)
		return functypeFLOOR;
	else if(name.compare("ceil")==0||name.compare("CEIL")==0)
		return functypeCEIL;
	//else if(name.compare("trunc")==0||name.compare("TRUNC")==0)
	//	return functypeTRUNC;
	else if(name.compare("round")==0||name.compare("ROUND")==0)
		return functypeROUND;
	else if(name.compare("sqrt")==0||name.compare("SQRT")==0)
		return functypeSQRT;
	else if(name.compare("ln")==0||name.compare("LN")==0
	     || name.compare("logE")==0||name.compare("LOGE")==0)
		return functypeLN;
	else if(name.compare("lb")==0||name.compare("LB")==0
	     || name.compare("log2")==0||name.compare("LOG2")==0)
		return functypeLB;
	else if(name.compare("lg")==0||name.compare("LG")==0
	     || name.compare("log10")==0||name.compare("LOG10")==0)
		return functypeLG;
	else if(name.compare("exp")==0||name.compare("EXP")==0)
		return functypeEXP;
	else if(name.compare("sin")==0||name.compare("SIN")==0)
		return functypeSIN;
	else if(name.compare("cos")==0||name.compare("COS")==0)
		return functypeCOS;
	else if(name.compare("tan")==0||name.compare("TAN")==0)
		return functypeTAN;
	else if(name.compare("asin")==0||name.compare("ASIN")==0)
		return functypeASIN;
	else if(name.compare("acos")==0||name.compare("ACOS")==0)
		return functypeACOS;
	else if(name.compare("atan")==0||name.compare("ATAN")==0)
		return functypeATAN;
	//
	else if(name.compare("min")==0||name.compare("MIN")==0)
		return functypeMIN;
	else if(name.compare("max")==0||name.compare("MAX")==0)
		return functypeMAX;
	else if(name.compare("mod")==0||name.compare("MOD")==0)
		return functypeMOD;
	else if(name.compare("pow")==0||name.compare("POW")==0)
		return functypePOW;
	//
	else if(name.compare("rcl")==0||name.compare("RCL")==0)
		return functypeRCL;
	//
	else if(name.compare("sto")==0||name.compare("STO")==0)
		return functypeSTO;
	//
	else
		return functypeNONE;
}

Expression::Tvaltype Expression::NodeFunction::valueType() const
{
	if(_functype==functypeNONE) return valtypeNONE;
	if(_children[1])
	{
		switch(_functype)
		{	// functions with 1 argument
			case functypeABS:
				return _children[1]->valueType();
			case functypeINT:
			case functypeFLOOR:
			case functypeCEIL:
			//case functypeTRUNC:
			case functypeROUND:
				return valtypeINT;
			case functypeFLOAT:
			case functypeSQRT:
			case functypeLN:
			case functypeLB:
			case functypeLG:
			case functypeEXP:
			case functypeSIN:
			case functypeCOS:
			case functypeTAN:
			case functypeASIN:
			case functypeACOS:
			case functypeATAN:
				return valtypeFLOAT;
			default: break;	// avoid compiler warnings of values not handled in switch
		}
		if(_children[0])
		{
			switch(_functype)
			{	// functions with 2 arguments
				case functypeMIN:
					if(_children[0]->evaluate()<_children[1]->evaluate())
						return _children[0]->valueType();
					else
						return _children[1]->valueType();
				case functypeMAX:
					if(_children[0]->evaluate()>_children[1]->evaluate())
						return _children[0]->valueType();
					else
						return _children[1]->valueType();
				case functypeMOD:
					if(_children[0]->valueType()==valtypeINT && _children[1]->valueType()==valtypeINT)
						return valtypeINT;
					else
						return valtypeFLOAT;
				case functypePOW:
					if(_children[0]->valueType()==valtypeINT && _children[1]->valueType()==valtypeINT)
						return valtypeINT;
					else
						return valtypeFLOAT;
				default: break;	// avoid compiler warnings of values not handled in switch
			}
		}
	}
	return valtypeNONE;
}

Expression::Value Expression::NodeFunction::evaluate() const
{
	if(_children[1])
	{
		switch(_functype)
		{	// functions with 1 argument
			case functypeNONE: return Value();
			case functypeABS:
				if(_children[1]->isInt())
					return Value(labs(Tint(_children[1]->evaluate())));
				else
					return Value(fabs(Tfloat(_children[1]->evaluate())));
			case functypeFLOAT: return Value(Tfloat(_children[1]->evaluate()));
			case functypeINT: return Value(Tint(Tfloat(_children[1]->evaluate())));
			case functypeFLOOR: return Value(Tint(floor(Tfloat(_children[1]->evaluate()))));
			case functypeCEIL:  return Value(Tint(ceil(Tfloat(_children[1]->evaluate()))));
			//case functypeTRUNC:  return Value(Tint(trunc(Tfloat(_children[1]->evaluate()))));
			case functypeROUND:  return Value(Tint(round(Tfloat(_children[1]->evaluate()))));
			case functypeSQRT:  return Value(sqrt(Tfloat(_children[1]->evaluate())));
			case functypeLN:   return Value(log(Tfloat(_children[1]->evaluate())));
			case functypeLB:  return Value(log2(Tfloat(_children[1]->evaluate())));
			case functypeLG: return Value(log10(Tfloat(_children[1]->evaluate())));
			case functypeEXP:   return Value(exp(Tfloat(_children[1]->evaluate())));
			case functypeSIN:   return Value(sin(Tfloat(_children[1]->evaluate())));
			case functypeCOS:   return Value(cos(Tfloat(_children[1]->evaluate())));
			case functypeTAN:   return Value(tan(Tfloat(_children[1]->evaluate())));
			case functypeASIN:   return Value(asin(Tfloat(_children[1]->evaluate())));
			case functypeACOS:   return Value(acos(Tfloat(_children[1]->evaluate())));
			case functypeATAN:   return Value(atan(Tfloat(_children[1]->evaluate())));
			default: break;	// avoid compiler warnings of values not handled in switch
		}
		if(_children[0])
		{
			switch(_functype)
			{	// functions with 2 arguments
				case functypeMIN:
					if(_children[0]->evaluate()<_children[1]->evaluate())
						return _children[0]->evaluate();
					else
						return _children[1]->evaluate();
				case functypeMAX:
					if(_children[0]->evaluate()>_children[1]->evaluate())
						return _children[0]->evaluate();
					else
						return _children[1]->evaluate();
				case functypeMOD:
					if(_children[0]->valueType()==valtypeINT && _children[1]->valueType()==valtypeINT)
						return Tint(_children[0]->evaluate()) % Tint(_children[1]->evaluate());
					else
						return Value(fmod(Tfloat(_children[0]->evaluate()),Tfloat(_children[1]->evaluate())));
				case functypePOW:
					if(_children[0]->valueType()==valtypeINT && _children[1]->valueType()==valtypeINT)
					{
						Tint b=_children[0]->evaluate();
						Tint e=_children[1]->evaluate();
						Tint p=1;
						for(Tint i=0;i<e;++i) p*=b;
						return Value(p);
					} else
						return Value(pow(Tfloat(_children[0]->evaluate()),Tfloat(_children[1]->evaluate())));
				default: break;	// avoid compiler warnings of values not handled in switch
			}
		}
	}
	return Value();
}

void Expression::NodeFunction::write(std::ostream& ostrm) const {
	switch(_functype)
	{	//
		case functypeNONE:	ostrm << "undef"; break;
		case functypeRCL:	ostrm << "RCL"; break;
		case functypeABS:	ostrm << "ABS"; break;
		case functypeFLOAT:	ostrm << "FLOAT"; break;
		case functypeINT:	ostrm << "INT"; break;
		case functypeFLOOR:	ostrm << "FLOOR"; break;
		case functypeCEIL:	ostrm << "CEIL"; break;
		//case functypeTRUNC:	ostrm << "TRUNC"; break;
		case functypeROUND:	ostrm << "ROUND"; break;
		case functypeSQRT:	ostrm << "SQRT"; break;
		case functypeLN:	ostrm << "LN"; break;
		case functypeLB:	ostrm << "LB"; break;
		case functypeLG:	ostrm << "LG"; break;
		case functypeEXP:	ostrm << "EXP"; break;
		case functypeSIN:	ostrm << "SIN"; break;
		case functypeCOS:	ostrm << "COS"; break;
		case functypeTAN:	ostrm << "TAN"; break;
		case functypeASIN:	ostrm << "ASIN"; break;
		case functypeACOS:	ostrm << "ACOS"; break;
		case functypeATAN:	ostrm << "ATAN"; break;
		//
		case functypeSTO:	ostrm << "STO"; break;
		case functypeMIN:	ostrm << "MIN"; break;
		case functypeMAX:	ostrm << "MAX"; break;
		case functypeMOD:	ostrm << "MOD"; break;
		case functypePOW:	ostrm << "POW"; break;
		default: break;	// avoid compiler warnings of values not handled in switch
	}
}


Expression::Tvaltype Expression::NodeFunctionVarSet::valueType() const
{
	if(_functype==functypeNONE) return valtypeNONE;
	if(_children[1])
	{
		switch(_functype)
		{	// functions with 1 argument
			case functypeRCL:
			{
				std::string varname("_localstore:"+static_cast<std::string>(*_children[1]));
				Expression::Tvaltype valtype = valtypeNONE;
				Expression::Variable* var=NULL;
				if(_variableset) var=_variableset->getVariable(varname);
				if(var) valtype=var->getType();
				return valtype;
			}
			default: break;	// avoid compiler warnings of values not handled in switch
		}
		if(_children[0])
		{
			switch(_functype)
			{	// functions with 2 arguments
				case functypeSTO:
					return _children[0]->valueType();
				default: break;	// avoid compiler warnings of values not handled in switch
			}
		}
	}
	return valtypeNONE;
}

Expression::Value Expression::NodeFunctionVarSet::evaluate() const
{
	if(_children[1])
	{
		switch(_functype)
		{	// functions with 1 argument
			case functypeRCL:
			{
				std::string varname("_localstore:"+static_cast<std::string>(*_children[1]));
				Expression::Value val;
				Expression::Variable* var=NULL;
				if(_variableset) var=_variableset->getVariable(varname);
				if(var) val=var->getValue();
				return val;
			}
			default: break;	// avoid compiler warnings of values not handled in switch
		}
		if(_children[0])
		{
			switch(_functype)
			{	// functions with 2 arguments
				case functypeSTO:
				{
					Expression::Value val=_children[0]->evaluate();
					std::string varname("_localstore:"+static_cast<std::string>(*_children[1]));
					if(_variableset) _variableset->setVariable(varname,val);
					return val;
				}
				default: break;	// avoid compiler warnings of values not handled in switch
			}
		}
	}
	return Value();
}

void Expression::NodeFunctionVarSet::write(std::ostream& ostrm) const {
	switch(_functype)
	{	//
		case functypeRCL:	ostrm << "RCL"; break;
		//
		case functypeSTO:	ostrm << "STO"; break;
		default: break;	// avoid compiler warnings of values not handled in switch
	}
}



void Expression::initializeRPN(const std::string& exprstr, bool genlabel)
{
	clear();
	std::stack<Node*> nodestack;
	// split expr to generate Nodes
	size_t startpos=0,endpos;
	while(startpos<exprstr.size()) {
		endpos=exprstr.find_first_of(" \t",startpos);
		if(endpos==std::string::npos)
		{
			endpos=exprstr.size()-1;
		} else if(endpos==startpos) {
			startpos=endpos+1;
			continue;
		} else {
			--endpos;
		}
		std::string token=exprstr.substr(startpos,endpos-startpos+1);
		size_t colonpos=token.find(":");
		if(token.size()==1 && token.find_first_of("+-*/")==0)
		{	// operator ................................................
			char op=token[0];
			if(nodestack.size()>=2)
			{
				Node* node1=nodestack.top(); nodestack.pop();
				Node* node2=nodestack.top(); nodestack.pop();
				Node* node0=new NodeOperation2(op,node2,node1);
				nodestack.push(node0);
				node1->setParent(node0);
				node2->setParent(node0);
			}
		} else if(token.find_first_not_of("0123456789.-E")==std::string::npos){
			// constant ................................................
			std::istringstream iss(token);
			if(token.find_first_not_of("0123456789-")==std::string::npos)
			{ // Tint
				Tint valInt=0;
				iss>>valInt;
				nodestack.push(new NodeConstant(valInt));
			}
			else
			{ // Tfloat
				Tfloat valFloat=0.;
				iss>>valFloat;
				nodestack.push(new NodeConstant(valFloat));
			}
			iss.str(std::string());
		} else if(colonpos!=std::string::npos) {
			// variable ................................................
			nodestack.push(new NodeVariable(_variableset->addVariable(token)));
		} else {
			// function ................................................
			//string funcname=token;
			NodeFunction::Efunctype func=NodeFunction::functype(token);
			if(func>NodeFunction::functypeMarkerVarSet2Arg && nodestack.size()>=2)
			{
				Node* node1=nodestack.top(); nodestack.pop();
				Node* node2=nodestack.top(); nodestack.pop();
				Node* node0=new NodeFunctionVarSet(func,_variableset,node1,node2);
				nodestack.push(node0);
				node1->setParent(node0);
			}
			else if(func>NodeFunction::functypeMarkerVarSet1Arg && nodestack.size()>=1)
			{
				Node* node1=nodestack.top(); nodestack.pop();
				Node* node0=new NodeFunctionVarSet(func,_variableset,node1);
				nodestack.push(node0);
				node1->setParent(node0);
			}
			else if(func>NodeFunction::functypeMarker2Arg && nodestack.size()>=2)
			{
				Node* node1=nodestack.top(); nodestack.pop();
				Node* node2=nodestack.top(); nodestack.pop();
				Node* node0=new NodeFunction(func,node1,node2);
				nodestack.push(node0);
				node1->setParent(node0);
			}
			else if(func>NodeFunction::functypeMarker1Arg && nodestack.size()>=1)
			{
				Node* node1=nodestack.top(); nodestack.pop();
				Node* node0=new NodeFunction(func,node1);
				nodestack.push(node0);
				node1->setParent(node0);
			}
		}
		startpos=endpos+2;
	}
	if(!nodestack.empty()) _rootnode=nodestack.top();
	if(genlabel) genLabel();
}

