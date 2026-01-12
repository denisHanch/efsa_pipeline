"""
Unit tests for the settings module.

Tests BaseSettings class including:
- Immutable update pattern
- Deep copy functionality
- Dictionary conversion (to_dict/from_dict)
- Validation of setting names
- String representations
"""

import pytest
from dataclasses import dataclass
from typing import Optional

from validation_pkg.utils.base_settings import BaseSettings


# ===== Test Fixtures =====

@dataclass
class SimpleSettings(BaseSettings):
    """Simple test settings class with basic fields."""
    param1: str = "default1"
    param2: int = 42
    param3: bool = True


@dataclass
class ComplexSettings(BaseSettings):
    """Complex test settings class with various types."""
    string_param: str = "test"
    int_param: int = 100
    float_param: float = 3.14
    bool_param: bool = False
    optional_param: Optional[str] = None
    list_param: list = None

    def __post_init__(self):
        if self.list_param is None:
            self.list_param = []


# ===== Tests =====

class TestBasicFunctionality:
    """Test basic BaseSettings functionality."""

    def test_create_settings(self):
        """Test creating settings with default values."""
        settings = SimpleSettings()
        assert settings.param1 == "default1"
        assert settings.param2 == 42
        assert settings.param3 is True

    def test_create_settings_with_custom_values(self):
        """Test creating settings with custom values."""
        settings = SimpleSettings(param1="custom", param2=100, param3=False)
        assert settings.param1 == "custom"
        assert settings.param2 == 100
        assert settings.param3 is False

    def test_settings_are_dataclasses(self):
        """Test that settings behave as dataclasses."""
        settings1 = SimpleSettings()
        settings2 = SimpleSettings()
        assert settings1 == settings2

        settings3 = SimpleSettings(param1="different")
        assert settings1 != settings3


class TestCopyMethod:
    """Test the copy() method."""

    def test_copy_creates_new_instance(self):
        """Test that copy creates a new instance."""
        original = SimpleSettings(param1="original")
        copied = original.copy()

        assert original is not copied
        assert original == copied

    def test_copy_is_deep_copy(self):
        """Test that copy is a deep copy."""
        original = ComplexSettings(list_param=[1, 2, 3])
        copied = original.copy()

        # Modify the copied list
        copied.list_param.append(4)

        # Original should be unchanged
        assert original.list_param == [1, 2, 3]
        assert copied.list_param == [1, 2, 3, 4]

    def test_copy_preserves_all_fields(self):
        """Test that copy preserves all field values."""
        original = ComplexSettings(
            string_param="test",
            int_param=999,
            float_param=2.71,
            bool_param=True,
            optional_param="optional"
        )
        copied = original.copy()

        assert copied.string_param == "test"
        assert copied.int_param == 999
        assert copied.float_param == 2.71
        assert copied.bool_param is True
        assert copied.optional_param == "optional"


class TestUpdateMethod:
    """Test the update() method (immutable pattern)."""

    def test_update_returns_new_instance(self):
        """Test that update returns a new instance."""
        original = SimpleSettings()
        updated = original.update(param1="new_value")

        assert original is not updated
        assert original.param1 == "default1"
        assert updated.param1 == "new_value"

    def test_update_single_parameter(self):
        """Test updating a single parameter."""
        settings = SimpleSettings()
        new_settings = settings.update(param1="updated")

        assert new_settings.param1 == "updated"
        assert new_settings.param2 == 42  # unchanged
        assert new_settings.param3 is True  # unchanged

    def test_update_multiple_parameters(self):
        """Test updating multiple parameters at once."""
        settings = SimpleSettings()
        new_settings = settings.update(param1="new", param2=999, param3=False)

        assert new_settings.param1 == "new"
        assert new_settings.param2 == 999
        assert new_settings.param3 is False

    def test_update_with_no_parameters(self):
        """Test update with no parameters returns copy."""
        original = SimpleSettings()
        updated = original.update()

        assert original == updated
        assert original is not updated

    def test_update_preserves_unchanged_fields(self):
        """Test that update preserves fields that weren't changed."""
        original = ComplexSettings(
            string_param="original",
            int_param=100,
            bool_param=True
        )
        updated = original.update(string_param="updated")

        assert updated.string_param == "updated"
        assert updated.int_param == 100  # preserved
        assert updated.bool_param is True  # preserved

    def test_update_raises_error_for_unknown_field(self):
        """Test that update raises ValueError for unknown fields."""
        settings = SimpleSettings()

        with pytest.raises(ValueError) as exc_info:
            settings.update(unknown_param="value")

        assert "Unknown setting(s) 'unknown_param'" in str(exc_info.value)
        assert "SimpleSettings" in str(exc_info.value)

    def test_update_raises_error_for_multiple_unknown_fields(self):
        """Test error message for multiple unknown fields."""
        settings = SimpleSettings()

        with pytest.raises(ValueError) as exc_info:
            settings.update(unknown1="val1", unknown2="val2", param1="valid")

        error_msg = str(exc_info.value)
        assert "unknown1" in error_msg
        assert "unknown2" in error_msg
        assert "Unknown setting(s)" in error_msg

    def test_update_error_lists_allowed_settings(self):
        """Test that error message lists allowed settings."""
        settings = SimpleSettings()

        with pytest.raises(ValueError) as exc_info:
            settings.update(invalid="value")

        error_msg = str(exc_info.value)
        assert "Allowed settings:" in error_msg
        assert "param1" in error_msg
        assert "param2" in error_msg
        assert "param3" in error_msg

    def test_update_chain(self):
        """Test chaining multiple update calls."""
        settings = SimpleSettings()
        result = (settings
                  .update(param1="step1")
                  .update(param2=100)
                  .update(param3=False))

        assert result.param1 == "step1"
        assert result.param2 == 100
        assert result.param3 is False
        assert settings.param1 == "default1"  # original unchanged


class TestToDictMethod:
    """Test the to_dict() method."""

    def test_to_dict_returns_dictionary(self):
        """Test that to_dict returns a dictionary."""
        settings = SimpleSettings()
        result = settings.to_dict()

        assert isinstance(result, dict)

    def test_to_dict_contains_all_fields(self):
        """Test that to_dict contains all fields."""
        settings = SimpleSettings(param1="test", param2=100, param3=False)
        result = settings.to_dict()

        assert result == {
            'param1': 'test',
            'param2': 100,
            'param3': False
        }

    def test_to_dict_with_complex_settings(self):
        """Test to_dict with complex settings."""
        settings = ComplexSettings(
            string_param="test",
            int_param=42,
            float_param=1.23,
            bool_param=True,
            optional_param="value",
            list_param=[1, 2, 3]
        )
        result = settings.to_dict()

        assert result == {
            'string_param': 'test',
            'int_param': 42,
            'float_param': 1.23,
            'bool_param': True,
            'optional_param': 'value',
            'list_param': [1, 2, 3]
        }

    def test_to_dict_with_none_values(self):
        """Test to_dict preserves None values."""
        settings = ComplexSettings(optional_param=None)
        result = settings.to_dict()

        assert 'optional_param' in result
        assert result['optional_param'] is None


class TestFromDictMethod:
    """Test the from_dict() class method."""

    def test_from_dict_creates_instance(self):
        """Test that from_dict creates a settings instance."""
        data = {'param1': 'test', 'param2': 99, 'param3': False}
        settings = SimpleSettings.from_dict(data)

        assert isinstance(settings, SimpleSettings)
        assert settings.param1 == 'test'
        assert settings.param2 == 99
        assert settings.param3 is False

    def test_from_dict_with_partial_data(self):
        """Test from_dict with partial data uses defaults."""
        data = {'param1': 'custom'}
        settings = SimpleSettings.from_dict(data)

        assert settings.param1 == 'custom'
        assert settings.param2 == 42  # default
        assert settings.param3 is True  # default

    def test_from_dict_with_empty_dict(self):
        """Test from_dict with empty dict uses all defaults."""
        settings = SimpleSettings.from_dict({})

        assert settings.param1 == 'default1'
        assert settings.param2 == 42
        assert settings.param3 is True

    def test_from_dict_raises_error_for_unknown_fields(self):
        """Test that from_dict raises ValueError for unknown fields."""
        data = {'param1': 'test', 'unknown_param': 'value'}

        with pytest.raises(ValueError) as exc_info:
            SimpleSettings.from_dict(data)

        assert "Unknown setting(s) 'unknown_param'" in str(exc_info.value)

    def test_from_dict_error_lists_allowed_settings(self):
        """Test that error message lists allowed settings."""
        data = {'invalid': 'value'}

        with pytest.raises(ValueError) as exc_info:
            SimpleSettings.from_dict(data)

        error_msg = str(exc_info.value)
        assert "Allowed settings:" in error_msg
        assert "param1" in error_msg

    def test_roundtrip_to_dict_from_dict(self):
        """Test roundtrip conversion to dict and back."""
        original = ComplexSettings(
            string_param="roundtrip",
            int_param=777,
            float_param=9.99,
            bool_param=True,
            optional_param="test"
        )

        # Convert to dict and back
        data = original.to_dict()
        restored = ComplexSettings.from_dict(data)

        assert original == restored


class TestStringRepresentations:
    """Test __str__ and __repr__ methods."""

    def test_str_format(self):
        """Test __str__ produces readable output."""
        settings = SimpleSettings(param1="test", param2=100, param3=False)
        result = str(settings)

        assert "SimpleSettings:" in result
        assert "param1: test" in result
        assert "param2: 100" in result
        assert "param3: False" in result

    def test_str_multiline(self):
        """Test that __str__ produces multiline output."""
        settings = SimpleSettings()
        result = str(settings)

        lines = result.split('\n')
        assert len(lines) > 1
        assert lines[0] == "SimpleSettings:"

    def test_repr_format(self):
        """Test __repr__ produces valid Python representation."""
        settings = SimpleSettings(param1="test", param2=100, param3=False)
        result = repr(settings)

        assert "SimpleSettings(" in result
        assert "param1='test'" in result
        assert "param2=100" in result
        assert "param3=False" in result

    def test_repr_can_be_evaluated(self):
        """Test that repr output looks like valid Python."""
        settings = SimpleSettings(param1="test", param2=100, param3=False)
        result = repr(settings)

        # Should start with class name and opening paren
        assert result.startswith("SimpleSettings(")
        assert result.endswith(")")


class TestImmutabilityPattern:
    """Test that the immutable pattern works correctly."""

    def test_original_unchanged_after_update(self):
        """Test that original settings are never modified."""
        original = SimpleSettings(param1="original", param2=1, param3=True)

        # Make multiple updates
        _ = original.update(param1="updated1")
        _ = original.update(param2=2)
        _ = original.update(param3=False)

        # Original should be completely unchanged
        assert original.param1 == "original"
        assert original.param2 == 1
        assert original.param3 is True

    def test_multiple_branches_from_same_original(self):
        """Test creating multiple branches from the same original."""
        base = SimpleSettings(param1="base", param2=0, param3=True)

        branch1 = base.update(param1="branch1")
        branch2 = base.update(param2=100)
        branch3 = base.update(param3=False)

        # Each branch should have only its changes
        assert branch1.param1 == "branch1"
        assert branch1.param2 == 0
        assert branch1.param3 is True

        assert branch2.param1 == "base"
        assert branch2.param2 == 100
        assert branch2.param3 is True

        assert branch3.param1 == "base"
        assert branch3.param2 == 0
        assert branch3.param3 is False

        # Base should be unchanged
        assert base.param1 == "base"
        assert base.param2 == 0
        assert base.param3 is True

    def test_shared_mutable_objects_are_deep_copied(self):
        """Test that mutable objects are deep copied."""
        original = ComplexSettings(list_param=[1, 2, 3])
        updated = original.update(string_param="updated")

        # Modify the list in updated
        updated.list_param.append(4)

        # Original list should be unchanged
        assert original.list_param == [1, 2, 3]
        assert updated.list_param == [1, 2, 3, 4]


class TestEdgeCases:
    """Test edge cases and special scenarios."""

    def test_settings_with_all_none_values(self):
        """Test settings where all optional values are None."""
        settings = ComplexSettings(
            string_param="test",
            int_param=0,
            float_param=0.0,
            bool_param=False,
            optional_param=None,
            list_param=None
        )

        # Should be able to convert to dict and back
        data = settings.to_dict()
        restored = ComplexSettings.from_dict(data)
        assert settings == restored

    def test_settings_equality(self):
        """Test that settings with same values are equal."""
        settings1 = SimpleSettings(param1="test", param2=100, param3=True)
        settings2 = SimpleSettings(param1="test", param2=100, param3=True)

        assert settings1 == settings2
        assert settings1 is not settings2

    def test_settings_inequality(self):
        """Test that settings with different values are not equal."""
        settings1 = SimpleSettings(param1="test1")
        settings2 = SimpleSettings(param1="test2")

        assert settings1 != settings2

    def test_update_with_same_values_creates_equal_object(self):
        """Test that updating with same values creates equal object."""
        original = SimpleSettings(param1="test", param2=100)
        updated = original.update(param1="test", param2=100)

        assert original == updated
        assert original is not updated  # But they're still different objects

    def test_empty_to_dict(self):
        """Test to_dict on settings with defaults."""
        settings = SimpleSettings()
        result = settings.to_dict()

        assert len(result) == 3
        assert all(key in result for key in ['param1', 'param2', 'param3'])


class TestValidatorIntegration:
    """Test BaseSettings with actual validator settings patterns."""

    @dataclass
    class MockValidatorSettings(BaseSettings):
        """Mock validator settings similar to actual validators."""
        validation_level: str = 'strict'
        coding_type: Optional[str] = None
        output_filename_suffix: Optional[str] = None
        output_subdir_name: Optional[str] = None
        check_duplicates: bool = True
        min_length: int = 100

    def test_validator_settings_creation(self):
        """Test creating settings like validators do."""
        settings = self.MockValidatorSettings()

        assert settings.validation_level == 'strict'
        assert settings.coding_type is None
        assert settings.check_duplicates is True
        assert settings.min_length == 100

    def test_validator_settings_update_pattern(self):
        """Test the common validator update pattern."""
        settings = self.MockValidatorSettings()
        settings = settings.update(
            validation_level='trust',
            coding_type='gz',
            output_filename_suffix='validated',
            min_length=500
        )

        assert settings.validation_level == 'trust'
        assert settings.coding_type == 'gz'
        assert settings.output_filename_suffix == 'validated'
        assert settings.min_length == 500

    def test_validator_settings_from_config(self):
        """Test loading settings from config dict."""
        config = {
            'validation_level': 'minimal',
            'coding_type': 'bz2',
            'check_duplicates': False
        }

        settings = self.MockValidatorSettings.from_dict(config)

        assert settings.validation_level == 'minimal'
        assert settings.coding_type == 'bz2'
        assert settings.check_duplicates is False
        assert settings.min_length == 100  # default


class TestTypeHints:
    """Test that type hints are preserved."""

    def test_fields_have_correct_types(self):
        """Test that field types match their annotations."""
        from dataclasses import fields

        simple_fields = {f.name: f.type for f in fields(SimpleSettings)}
        assert simple_fields['param1'] == str
        assert simple_fields['param2'] == int
        assert simple_fields['param3'] == bool

    def test_optional_fields(self):
        """Test that optional fields work correctly."""
        from dataclasses import fields

        complex_fields = {f.name: f.type for f in fields(ComplexSettings)}
        assert 'Optional' in str(complex_fields['optional_param'])
