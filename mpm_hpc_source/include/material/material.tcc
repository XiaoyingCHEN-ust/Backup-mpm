//! Get material property
template <unsigned Tdim>
template <typename Ttype>
Ttype mpm::Material<Tdim>::property(const std::string& key) {
  try {
    return properties_[key].template get<Ttype>();
  } catch (std::exception& except) {
    console_->error("Property call to material parameter not found: {}",
                    except.what());
    throw std::runtime_error(
        "Property call to material parameter not found or invalid type");
  }
}

//! Get optional material property
template <unsigned Tdim>
template <typename Ttype>
Ttype mpm::Material<Tdim>::property_or(const std::string& key,
                                       const Ttype& default_value) {
  try {
    auto itr = properties_.find(key);
    if (itr == properties_.end() || itr->is_null()) return default_value;
    return itr->template get<Ttype>();
  } catch (std::exception& except) {
    console_->error("Optional property call to material parameter invalid: {}",
                    except.what());
    throw std::runtime_error("Optional material parameter invalid type");
  }
}
