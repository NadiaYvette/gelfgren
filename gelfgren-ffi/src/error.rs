//! Error handling for FFI layer.
//!
//! Errors are stored in thread-local storage and can be retrieved via
//! `gelfgren_last_error_message()`.

use std::cell::RefCell;
use std::ffi::CString;
use std::os::raw::c_char;
use std::panic;

use gelfgren_core::error::GelfgrenError;
use crate::types::GelfgrenErrorCode;

thread_local! {
    static LAST_ERROR: RefCell<Option<String>> = RefCell::new(None);
}

/// Sets the last error message for the current thread.
pub(crate) fn set_last_error(error: String) {
    LAST_ERROR.with(|e| {
        *e.borrow_mut() = Some(error);
    });
}

/// Clears the last error message for the current thread.
pub(crate) fn clear_last_error() {
    LAST_ERROR.with(|e| {
        *e.borrow_mut() = None;
    });
}

/// Retrieves the last error message for the current thread.
///
/// Returns NULL if no error has occurred.
///
/// # Safety
///
/// The returned string is valid until the next call to any Gelfgren function
/// on this thread. The caller must NOT free the returned pointer.
///
/// # Example (C)
///
/// ```c
/// const char* error = gelfgren_last_error_message();
/// if (error != NULL) {
///     fprintf(stderr, "Error: %s\n", error);
/// }
/// ```
#[no_mangle]
pub unsafe extern "C" fn gelfgren_last_error_message() -> *const c_char {
    LAST_ERROR.with(|e| {
        match &*e.borrow() {
            Some(err) => {
                // Allocate a new CString (will leak, but that's acceptable for error reporting)
                match CString::new(err.as_str()) {
                    Ok(c_str) => c_str.into_raw() as *const c_char,
                    Err(_) => std::ptr::null(),
                }
            }
            None => std::ptr::null(),
        }
    })
}

/// Clears the last error message.
///
/// # Example (C)
///
/// ```c
/// gelfgren_clear_last_error();
/// ```
#[no_mangle]
pub extern "C" fn gelfgren_clear_last_error() {
    clear_last_error();
}

/// Converts a Rust Result to an error code, storing the error message on failure.
pub(crate) fn result_to_code<T>(result: Result<T, GelfgrenError>) -> (Option<T>, GelfgrenErrorCode) {
    match result {
        Ok(val) => {
            clear_last_error();
            (Some(val), GelfgrenErrorCode::Success)
        }
        Err(e) => {
            let (code, msg) = error_to_code_and_message(&e);
            set_last_error(msg);
            (None, code)
        }
    }
}

/// Converts GelfgrenError to error code and message.
fn error_to_code_and_message(error: &GelfgrenError) -> (GelfgrenErrorCode, String) {
    match error {
        GelfgrenError::InvalidArgument(msg) => {
            (GelfgrenErrorCode::InvalidArgument, msg.clone())
        }
        GelfgrenError::InvalidDimension(msg) => {
            (GelfgrenErrorCode::InvalidDimension, msg.clone())
        }
        GelfgrenError::DivisionByZero => {
            (GelfgrenErrorCode::DivisionByZero, "Division by zero".to_string())
        }
        GelfgrenError::SingularMatrix => {
            (GelfgrenErrorCode::SingularMatrix, "Singular matrix".to_string())
        }
        GelfgrenError::NotImplemented(msg) => {
            (GelfgrenErrorCode::NotImplemented, msg.clone())
        }
        GelfgrenError::ConvergenceFailure(msg) => {
            (GelfgrenErrorCode::Unknown, format!("Convergence failure: {}", msg))
        }
        GelfgrenError::Overflow => {
            (GelfgrenErrorCode::Unknown, "Numerical overflow".to_string())
        }
        GelfgrenError::Underflow => {
            (GelfgrenErrorCode::Unknown, "Numerical underflow".to_string())
        }
    }
}

/// Macro to catch panics and convert to FFI errors.
///
/// Usage:
/// ```rust,ignore
/// ffi_catch! {
///     // Code that might panic
///     some_operation()?;
/// }
/// ```
#[macro_export]
macro_rules! ffi_catch {
    ($($body:tt)*) => {
        match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| { $($body)* })) {
            Ok(result) => result,
            Err(_) => {
                $crate::error::set_last_error("Panic occurred in Rust code".to_string());
                $crate::types::GelfgrenErrorCode::Unknown
            }
        }
    };
}

/// Macro to check for null pointers and return error code if null.
#[macro_export]
macro_rules! check_null {
    ($ptr:expr, $name:expr) => {
        if $ptr.is_null() {
            $crate::error::set_last_error(format!("Null pointer: {}", $name));
            return $crate::types::GelfgrenErrorCode::NullPointer;
        }
    };
}

/// Macro to check for null pointers and return null pointer if null (for pointer-returning functions).
#[macro_export]
macro_rules! check_null_ptr {
    ($ptr:expr, $name:expr) => {
        if $ptr.is_null() {
            $crate::error::set_last_error(format!("Null pointer: {}", $name));
            return std::ptr::null_mut();
        }
    };
}

/// Macro to check array dimensions and return error code if invalid.
#[macro_export]
macro_rules! check_dimension {
    ($size:expr, $min:expr, $name:expr) => {
        if $size < $min {
            $crate::error::set_last_error(format!(
                "Invalid dimension for {}: {} (minimum: {})",
                $name, $size, $min
            ));
            return $crate::types::GelfgrenErrorCode::InvalidDimension;
        }
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_storage() {
        clear_last_error();
        set_last_error("Test error".to_string());

        unsafe {
            let msg = gelfgren_last_error_message();
            assert!(!msg.is_null());
            let c_str = std::ffi::CStr::from_ptr(msg);
            assert_eq!(c_str.to_str().unwrap(), "Test error");
        }
    }

    #[test]
    fn test_error_conversion() {
        let err = GelfgrenError::InvalidArgument("Bad input".to_string());
        let (code, msg) = error_to_code_and_message(&err);
        assert_eq!(code, GelfgrenErrorCode::InvalidArgument);
        assert_eq!(msg, "Bad input");
    }

    #[test]
    fn test_result_to_code_success() {
        let result: Result<i32, GelfgrenError> = Ok(42);
        let (val, code) = result_to_code(result);
        assert_eq!(val, Some(42));
        assert_eq!(code, GelfgrenErrorCode::Success);
    }

    #[test]
    fn test_result_to_code_error() {
        let result: Result<i32, GelfgrenError> = Err(GelfgrenError::DivisionByZero);
        let (val, code) = result_to_code(result);
        assert_eq!(val, None);
        assert_eq!(code, GelfgrenErrorCode::DivisionByZero);
    }
}
