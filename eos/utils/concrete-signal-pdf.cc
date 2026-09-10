/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015-2026 Danny van Dyk
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <eos/maths/integrate.hh>
#include <eos/utils/concrete-signal-pdf.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

namespace eos
{
    ConcreteSignalPDF::ConcreteSignalPDF(const QualifiedName & name, const Parameters & parameters, const Kinematics & kinematics, const Options & options,
                                         const QualifiedName & unnormalized_pdf, const QualifiedName & normalization) :
        _name(name),
        _parameters(parameters),
        _kinematics(kinematics),
        _options(options),
        _unnormalized_pdf(Observable::make(unnormalized_pdf, parameters, kinematics, options)),
        _normalization(Observable::make(normalization, parameters, kinematics, options)),
        _integration_variable(""),
        _padding_fraction(0.0)
    {
        if (_unnormalized_pdf == nullptr)
        {
            throw InternalError("ConcreteSignalPDF: failed to construct unnormalized pdf from " + unnormalized_pdf.str());
        }

        if (_normalization == nullptr)
        {
            throw InternalError("ConcreteSignalPDF: failed to construct normalization from " + normalization.str());
        }
    }

    ConcreteSignalPDF::ConcreteSignalPDF(const QualifiedName & name, const Parameters & parameters, const Kinematics & kinematics, const Options & options,
                                         const QualifiedName & unnormalized_pdf, const std::string & integration_variable, const double & padding_fraction) :
        _name(name),
        _parameters(parameters),
        _kinematics(kinematics),
        _options(options),
        _unnormalized_pdf(Observable::make(unnormalized_pdf, parameters, kinematics, options)),
        _normalization(nullptr),
        _integration_variable(integration_variable),
        _padding_fraction(padding_fraction)
    {
        if (_unnormalized_pdf == nullptr)
        {
            throw InternalError("ConcreteSignalPDF: failed to construct unnormalized pdf from " + unnormalized_pdf.str());
        }
    }

    const QualifiedName &
    ConcreteSignalPDF::name() const
    {
        return _name;
    }

    double
    ConcreteSignalPDF::evaluate() const
    {
        if (auto result = _unnormalized_pdf->evaluate(); result > 0.0) [[likely]]
        {
            return std::log(result);
        }

        return -std::numeric_limits<double>::infinity();
    }

    double
    ConcreteSignalPDF::evaluate_linear() const
    {
        if (auto result = _unnormalized_pdf->evaluate(); result > 0.0) [[likely]]
        {
            return result;
        }

        return 0.0;
    }

    double
    ConcreteSignalPDF::normalization() const
    {
        if (_normalization == nullptr)
        {
            // Numerical integration of the PDF over the kinematic range
            // corrected to be shrinked by the padding fraction

            const std::string variable_name = _integration_variable;

            const double v_min_padded = _kinematics[variable_name + "_min"].evaluate();
            const double v_max_padded = _kinematics[variable_name + "_max"].evaluate();

            const double v_min = (1 + _padding_fraction) / (1 + 2 * _padding_fraction) * v_min_padded + _padding_fraction / (1 + 2 * _padding_fraction) * v_max_padded;
            const double v_max = v_min * (1 + _padding_fraction) / _padding_fraction - v_min_padded / _padding_fraction;

            Kinematics kinematics = _kinematics.clone();
            auto       v          = kinematics[variable_name];

            auto unnormalized_function = Observable::make(_unnormalized_pdf->name(), _parameters, kinematics, _options);

            const auto integrand = [&](const double & x) -> double
            {
                v = x;
                return unnormalized_function->evaluate();
            };

            // const auto   cfg  = GSL::QAGS::Config().epsabs(0.0).epsrel(1e-4).key(1);
            // const double norm = integrate<GSL::QAGS>(integrand, v_min, v_max, cfg);

            const auto trapezoidal_integral = [&](const double a, const double b, const unsigned n_points) -> double
            {
                if (b <= a)
                {
                    return 0.0;
                }

                const double h = (b - a) / (n_points - 1);

                double integral = 0.0;

                for (unsigned i = 0; i < n_points; ++i)
                {
                    const double x = a + i * h;

                    v = x;

                    const double y = unnormalized_function->evaluate();

                    if (! std::isfinite(y))
                    {
                        std::cerr << "Bad PDF value at q2 = " << x << std::endl;
                        throw std::runtime_error("Non-finite PDF");
                    }

                    // trapezoidal weights
                    if (i == 0 || i == n_points - 1)
                    {
                        integral += 0.5 * y;
                    }
                    else
                    {
                        integral += y;
                    }
                }

                return integral * h;
            };

            double norm = 0.0;

            if (v_min < 9.0 && v_max > 10.1)
            {
                norm += trapezoidal_integral(v_min, 9.0, 125);
                norm += trapezoidal_integral(9.0, 10.1, 400);
                norm += trapezoidal_integral(10.1, v_max, 75);
            }

            if (v_max <= 9.0)
            {
                norm += trapezoidal_integral(v_min, v_max, 500);
            }

            if (norm > 0.0)
            {
                return std::log(norm);
            }

            return -std::numeric_limits<double>::infinity();
        }

        if (auto result = _normalization->evaluate(); result > 0.0) [[likely]]
        {
            return std::log(result);
        }

        return -std::numeric_limits<double>::infinity();
    }

    ObservablePtr
    ConcreteSignalPDF::unnormalized_pdf() const
    {
        return _unnormalized_pdf;
    }

    ObservablePtr
    ConcreteSignalPDF::normalization_observable() const
    {
        return _normalization;
    }

    Parameters
    ConcreteSignalPDF::parameters()
    {
        return _parameters;
    }

    Kinematics
    ConcreteSignalPDF::kinematics()
    {
        return _kinematics;
    }

    Options
    ConcreteSignalPDF::options()
    {
        return _options;
    }

    DensityPtr
    ConcreteSignalPDF::clone() const
    {
        if (_normalization == nullptr)
        {
            return DensityPtr(
                    new ConcreteSignalPDF(_name, _parameters.clone(), _kinematics.clone(), _options, _unnormalized_pdf->name(), _integration_variable, _padding_fraction));
        }
        return DensityPtr(new ConcreteSignalPDF(_name, _parameters.clone(), _kinematics.clone(), _options, _unnormalized_pdf->name(), _normalization->name()));
    }

    DensityPtr
    ConcreteSignalPDF::clone(const Parameters & parameters) const
    {
        if (_normalization == nullptr)
        {
            return DensityPtr(new ConcreteSignalPDF(_name, parameters, _kinematics.clone(), _options, _unnormalized_pdf->name(), _integration_variable, _padding_fraction));
        }
        return DensityPtr(new ConcreteSignalPDF(_name, parameters, _kinematics.clone(), _options, _unnormalized_pdf->name(), _normalization->name()));
    }

    Density::Iterator
    ConcreteSignalPDF::begin() const
    {
        return Density::Iterator(_descriptions.cbegin());
    }

    Density::Iterator
    ConcreteSignalPDF::end() const
    {
        return Density::Iterator(_descriptions.cend());
    }

    ConcreteSignalPDFEntry::ConcreteSignalPDFEntry(const QualifiedName & name, const std::string & description, const Options & default_options, const QualifiedName & numerator,
                                                   const QualifiedName & normalization, const std::vector<std::string> & numerator_kinematic_names,
                                                   const std::vector<std::string> & normalization_kinematic_names) :
        _name(name),
        _description(description),
        _default_options(default_options),
        _numerator(numerator),
        _normalization(normalization),
        _numerator_kinematic_names(numerator_kinematic_names),
        _normalization_kinematic_names(normalization_kinematic_names)
    {
    }

    ConcreteSignalPDFEntry::ConcreteSignalPDFEntry(const QualifiedName & name, const std::string & description, const Options & default_options, const QualifiedName & numerator,
                                                   const std::vector<std::string> & numerator_kinematic_names, const std::vector<std::string> & normalization_kinematic_names,
                                                   const double & padding_fraction) :
        _name(name),
        _description(description),
        _default_options(default_options),
        _numerator(numerator),
        _normalization(QualifiedName("null::null")),
        _numerator_kinematic_names(numerator_kinematic_names),
        _normalization_kinematic_names(normalization_kinematic_names),
        _padding_fraction(padding_fraction)
    {
    }

    ConcreteSignalPDFEntry::~ConcreteSignalPDFEntry() = default;

    const QualifiedName &
    ConcreteSignalPDFEntry::name() const
    {
        return _name;
    }

    const std::string &
    ConcreteSignalPDFEntry::description() const
    {
        return _description;
    }

    SignalPDFEntry::NumeratorKinematicVariableIterator
    ConcreteSignalPDFEntry::begin_numerator_kinematic_variables() const
    {
        return SignalPDFEntry::NumeratorKinematicVariableIterator(_numerator_kinematic_names.cbegin());
    }

    SignalPDFEntry::NumeratorKinematicVariableIterator
    ConcreteSignalPDFEntry::end_numerator_kinematic_variables() const
    {
        return SignalPDFEntry::NumeratorKinematicVariableIterator(_numerator_kinematic_names.cend());
    }

    SignalPDFEntry::DenominatorKinematicVariableIterator
    ConcreteSignalPDFEntry::begin_denominator_kinematic_variables() const
    {
        return SignalPDFEntry::DenominatorKinematicVariableIterator(_normalization_kinematic_names.cbegin());
    }

    SignalPDFEntry::DenominatorKinematicVariableIterator
    ConcreteSignalPDFEntry::end_denominator_kinematic_variables() const
    {
        return SignalPDFEntry::DenominatorKinematicVariableIterator(_normalization_kinematic_names.cend());
    }

    SignalPDFPtr
    ConcreteSignalPDFEntry::make(const Parameters & parameters, const Kinematics & kinematics, const Options & options) const
    {
        if (_normalization.str() == "null::null")
        {
            return SignalPDFPtr(
                    new ConcreteSignalPDF(_name, parameters, kinematics, _default_options + options, _numerator, _numerator_kinematic_names.front(), _padding_fraction));
        }
        return SignalPDFPtr(new ConcreteSignalPDF(_name, parameters, kinematics, _default_options + options, _numerator, _normalization));
    }

    std::ostream &
    ConcreteSignalPDFEntry::insert(std::ostream & os) const
    {
        return os;
    }

    template <> struct WrappedForwardIteratorTraits<SignalPDFEntry::NumeratorKinematicVariableIteratorTag>
    {
            using UnderlyingIterator = std::vector<std::string>::const_iterator;
    };
    template class WrappedForwardIterator<SignalPDFEntry::NumeratorKinematicVariableIteratorTag, const std::string &>;

    template <> struct WrappedForwardIteratorTraits<SignalPDFEntry::DenominatorKinematicVariableIteratorTag>
    {
            using UnderlyingIterator = std::vector<std::string>::const_iterator;
    };
    template class WrappedForwardIterator<SignalPDFEntry::DenominatorKinematicVariableIteratorTag, const std::string &>;
} // namespace eos
