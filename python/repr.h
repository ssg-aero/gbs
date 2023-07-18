#pragma once
#include <gbs-io/tojson.h>
using namespace gbs;

inline std::string build_rep(const auto &cls)
{
    rapidjson::StringBuffer buffer;
    rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
    rapidjson::Document d;
    d.SetObject();
    auto crv_val = make_json(cls, d.GetAllocator());
    crv_val.Accept(writer);
    return std::string(buffer.GetString());
}