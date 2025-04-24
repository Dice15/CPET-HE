#pragma once

#include "util/basis.h"
#include <cstdint>
#include <memory>


namespace cpet
{
	class TextParameter
	{
	public:
		TextParameter() :basis_(nullptr), next_param_(nullptr) {}

		std::shared_ptr<const Basis> basis() const;

		std::shared_ptr<const TextParameter> next_param() const;

		static void create_parameter_chain(
			Basis basis,
			std::shared_ptr<const TextParameter>& first_param,
			std::shared_ptr<const TextParameter>& last_param);

	private:
		static std::shared_ptr<const TextParameter> internal_create_parameter_chain(
			Basis& basis,
			std::shared_ptr<const TextParameter>& last_param);

		std::shared_ptr<const Basis> basis_;

		std::shared_ptr<const TextParameter> next_param_;
	};
}