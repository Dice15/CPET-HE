#pragma once

#include "context.h"
#include "secretkey.h"
#include "publickey.h"
#include "switchkey.h"
#include <memory>


namespace cpet
{
	class KeyGenerator
	{
	public:
		KeyGenerator();

		KeyGenerator(const Context& context);

		void secret_key(SecretKey& destination) const;

		void create_public_key(PublicKey& destination) const;

		void create_evaluate_key(SwitchKey& destination) const;

	private:
		void internal_create_switch_key(const cpet::RnsCycloRing& s1, const cpet::RnsCycloRing& s2, SwitchKey& destination) const;

		std::shared_ptr<Context> context_;

		SecretKey sk_;


	};
}