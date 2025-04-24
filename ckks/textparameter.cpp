
#include "textparameter.h"
#include <cstdint>
#include <memory>


namespace cpet
{
	std::shared_ptr<const Basis> TextParameter::basis() const
	{
		return basis_;
	}

	std::shared_ptr<const TextParameter> TextParameter::next_param() const
	{
		return next_param_;
	}

	void TextParameter::create_parameter_chain(
		Basis basis,
		std::shared_ptr<const TextParameter>& first_param,
		std::shared_ptr<const TextParameter>& last_param)
	{
		first_param = internal_create_parameter_chain(basis, last_param);
	}

	std::shared_ptr<const TextParameter> TextParameter::internal_create_parameter_chain(
		Basis& basis,
		std::shared_ptr<const TextParameter>& last_param)
	{
		if (basis.size() == 0)
		{
			return nullptr;
		}

		TextParameter param;

		param.basis_ = std::make_shared<const Basis>(basis);
		basis.drop_q();
		param.next_param_ = internal_create_parameter_chain(basis, last_param);

		if (param.next_param_ == nullptr)
		{
			return last_param = std::make_shared<const TextParameter>(param);
		}

		return std::make_shared<const TextParameter>(param);
	}
}