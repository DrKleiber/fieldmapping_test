//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FieldMappingValue.h"
// MOOSE includes
#include "MooseVariableFE.h"

// C++ includes
#include <numeric>

registerMooseObject("ExampleApp", FieldMappingValue);

InputParameters
FieldMappingValue::validParams()
{
  InputParameters params = ElementVariableVectorPostprocessor::validParams();

  params.addClassDescription("Samples values of elemental variable(s).");

  params += SamplerBase::validParams();

//  params.set<std::string>("built_by_action") = "add_user_object";

  return params;
}

FieldMappingValue::FieldMappingValue(const InputParameters & parameters)
  : ElementVariableVectorPostprocessor(parameters), SamplerBase(parameters, this, _communicator)
{
  // ensure that variables are elemental, i.e., not scalar and and not nodal
  for (unsigned int i = 0; i < _coupled_moose_vars.size(); i++)
    if (_coupled_moose_vars[i]->feType().family == SCALAR || _coupled_moose_vars[i]->isNodal())
      paramError(
          "variable", "The variable '", _coupled_moose_vars[i]->name(), "' is not elemental.");

  std::vector<std::string> var_names(_coupled_moose_vars.size());
  _values.resize(_coupled_moose_vars.size());

  for (unsigned int i = 0; i < _coupled_moose_vars.size(); i++)
    var_names[i] = _coupled_moose_vars[i]->name();

  // Initialize the data structures in SamplerBase
  SamplerBase::setupVariables(var_names);
}

const std::vector<Real> &
FieldMappingValue::getValue() const
{
  return _values;
}

void
FieldMappingValue::initialize()
{
  SamplerBase::initialize();
}

void
FieldMappingValue::execute()
{
  for (unsigned int i = 0; i < _coupled_moose_vars.size(); i++)
    _values[i] = _coupled_standard_moose_vars[i]->getElementalValue(_current_elem);

  SamplerBase::addSample(_current_elem->centroid(), _current_elem->id(), _values);
}

void
FieldMappingValue::finalize()
{
  SamplerBase::finalize();
}

void
FieldMappingValue::threadJoin(const UserObject & y)
{
  const FieldMappingValue & vpp = static_cast<const FieldMappingValue &>(y);

  SamplerBase::threadJoin(vpp);
}
