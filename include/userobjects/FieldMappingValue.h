//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementVariableVectorPostprocessor.h"
#include "SamplerBase.h"

#include "libmesh/mesh_tools.h"

/**
 * Computes the average value of a variable on each block
 */
class FieldMappingValue :  public ElementVariableVectorPostprocessor, protected SamplerBase
{
public:
  static InputParameters validParams();
  FieldMappingValue(const InputParameters & parameters);
  const std::vector<Real> & getValue() const;

  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;

  // Let the SamplerBase version of threadJoin() take part in the
  // overload resolution process, otherwise we get warnings about
  // overloaded virtual functions and "hiding" in debug mode.
  using SamplerBase::threadJoin;
  virtual void threadJoin(const UserObject & y) override;

protected:
  std::vector<Real> _values;

};
