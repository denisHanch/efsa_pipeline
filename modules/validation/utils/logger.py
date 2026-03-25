"""Logging configuration for the bioinformatics validation package."""

import sys
import re
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, List
from threading import Lock
from dataclasses import dataclass
import structlog
import logging
import time


@dataclass
class FileTimingSummary:
    """Simple timing summary for a single file."""
    input_file: str
    validator_type: str  # "genome", "read", "feature"
    elapsed_time: float


# ANSI color codes for console output
class Colors:
    """ANSI color codes for terminal output."""
    RESET = '\033[0m'
    BOLD = '\033[1m'
    RED = '\033[91m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    MAGENTA = '\033[95m'
    CYAN = '\033[96m'
    GRAY = '\033[90m'


def add_log_level_colors(_, level: str, event_dict: dict) -> dict:
    """Add colors to log level in console output."""
    level_colors = {
        "debug": Colors.GRAY,
        "info": Colors.BLUE,
        "warning": Colors.YELLOW,
        "error": Colors.RED,
        "critical": Colors.RED + Colors.BOLD,
    }

    level_lower = level.lower()
    color = level_colors.get(level_lower, "")
    if color:
        event_dict["level"] = f"{color}{level.upper()}{Colors.RESET}"
    else:
        event_dict["level"] = level.upper()

    return event_dict


def format_process_info(logger, method_name, event_dict: dict) -> dict:
    """Format process/worker information for display in console logs."""
    worker_id = event_dict.get("worker_id")
    file_context = event_dict.get("file_context")
    category = event_dict.get("category")

    context_parts = []

    if worker_id:
        context_parts.append(f"{Colors.CYAN}Worker-{worker_id}{Colors.RESET}")

    if file_context:
        if len(file_context) > 40:
            file_context = "..." + file_context[-37:]
        context_parts.append(f"{Colors.MAGENTA}{file_context}{Colors.RESET}")

    if category and category not in ['validation_pipeline']:
        category_color = {
            'genome': Colors.GREEN,
            'read': Colors.BLUE,
            'feature': Colors.YELLOW,
            'inter-file': Colors.CYAN
        }.get(category, Colors.GRAY)
        context_parts.append(f"{category_color}{category}{Colors.RESET}")

    if context_parts:
        event_dict["context"] = f"[{' '.join(context_parts)}]"

    return event_dict


def get_incremented_path(path: Path, separator: str = "_") -> Path:
    """Get next available filename by auto-incrementing if file exists."""
    path = Path(path)

    if not path.exists():
        return path

    stem = path.stem
    suffix = path.suffix
    parent = path.parent

    match = re.match(r'^(.+)_(\d+)$', stem)
    if match:
        base_stem = match.group(1)
        start_counter = int(match.group(2)) + 1
    else:
        base_stem = stem
        start_counter = 1

    counter = start_counter
    while True:
        new_name = f"{base_stem}{separator}{counter:03d}{suffix}"
        new_path = parent / new_name
        if not new_path.exists():
            return new_path
        counter += 1

        if counter > 9999:
            raise RuntimeError(f"Too many incremented files for {path}. Maximum is 9999.")


class ValidationLogger:
    """Logger for validation package with structured logging support."""

    def __init__(self):
        """Initialize logger."""
        self.log_file: Optional[Path] = None
        self.validation_issues = []
        self._issues_lock = Lock()
        # Eagerly grab the structlog logger so fallback instances work if
        # structlog has already been configured (e.g. by setup_logging in main.py).
        self.logger = structlog.get_logger("validation_pipeline")
        self._timers: Dict[str, float] = {}
        self._timers_lock = Lock()
        self.file_timings: List[FileTimingSummary] = []
        self._timings_lock = Lock()

    def setup(
        self,
        console_level: str = "INFO",
        log_file: Optional[Path] = None,
        clear_previous_issues: bool = True
    ):
        """Set up logging handlers."""
        if clear_previous_issues:
            self.clear_issues()
            self.clear_file_timings()

        if log_file:
            log_file = Path(log_file)
            log_file.parent.mkdir(parents=True, exist_ok=True)

            log_file = get_incremented_path(log_file)
            self.log_file = log_file

            processors = [
                structlog.contextvars.merge_contextvars,
                structlog.stdlib.add_log_level,
                structlog.stdlib.add_logger_name,
                structlog.processors.TimeStamper(fmt="iso"),
                structlog.processors.StackInfoRenderer(),
                structlog.processors.format_exc_info,
                structlog.stdlib.ProcessorFormatter.wrap_for_formatter,
            ]

            structlog.configure(
                processors=processors,
                wrapper_class=structlog.stdlib.BoundLogger,
                context_class=dict,
                logger_factory=structlog.stdlib.LoggerFactory(),
                cache_logger_on_first_use=False,
            )

            stdlib_logger = logging.getLogger("validation_pipeline")
            stdlib_logger.handlers.clear()
            stdlib_logger.setLevel(logging.DEBUG)
            stdlib_logger.propagate = False

            file_handler = logging.FileHandler(log_file, mode='w', encoding='utf-8')
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(
                structlog.stdlib.ProcessorFormatter(
                    processor=structlog.processors.JSONRenderer(),
                    foreign_pre_chain=processors,
                )
            )
            stdlib_logger.addHandler(file_handler)

            console_pre_chain = processors[:-1] + [
                format_process_info,
                add_log_level_colors,
                processors[-1],
            ]

            console_handler = logging.StreamHandler(sys.stdout)
            console_handler.setLevel(getattr(logging, console_level.upper()))
            console_handler.setFormatter(
                structlog.stdlib.ProcessorFormatter(
                    processor=structlog.dev.ConsoleRenderer(colors=True),
                    foreign_pre_chain=console_pre_chain,
                )
            )
            stdlib_logger.addHandler(console_handler)

        else:
            processors = [
                structlog.contextvars.merge_contextvars,
                structlog.processors.add_log_level,
                structlog.processors.TimeStamper(fmt="iso"),
                structlog.processors.StackInfoRenderer(),
                structlog.processors.format_exc_info,
            ]

            console_processors = processors + [
                format_process_info,
                add_log_level_colors,
                structlog.dev.ConsoleRenderer(colors=True)
            ]

            structlog.configure(
                processors=console_processors,
                wrapper_class=structlog.make_filtering_bound_logger(
                    getattr(structlog.stdlib.logging, console_level.upper(), structlog.stdlib.logging.INFO)
                ),
                context_class=dict,
                logger_factory=structlog.PrintLoggerFactory(file=sys.stdout),
                cache_logger_on_first_use=False,
            )

        self.logger = structlog.get_logger("validation_pipeline")

        if log_file:
            self.info(f"Detailed log file: {log_file}")

    def reconfigure_level(
        self,
        console_level: str = "INFO",
        enable_file_logging: bool = True,
        log_file: Optional[Path] = None
    ):
        """Reconfigure logging level after initial setup."""
        console_level_upper = console_level.upper()

        stdlib_logger = logging.getLogger("validation_pipeline")

        if stdlib_logger.handlers:
            for handler in stdlib_logger.handlers:
                if isinstance(handler, logging.StreamHandler) and not isinstance(handler, logging.FileHandler):
                    handler.setLevel(getattr(logging, console_level_upper))
                    self.debug(f"Updated console logging level to {console_level_upper}")

            if not enable_file_logging:
                stdlib_logger.handlers = [
                    h for h in stdlib_logger.handlers
                    if not isinstance(h, logging.FileHandler)
                ]
        else:
            processors = [
                structlog.contextvars.merge_contextvars,
                structlog.processors.add_log_level,
                structlog.processors.TimeStamper(fmt="iso"),
                structlog.processors.StackInfoRenderer(),
                structlog.processors.format_exc_info,
            ]

            console_processors = processors + [
                format_process_info,
                add_log_level_colors,
                structlog.dev.ConsoleRenderer(colors=True)
            ]

            structlog.configure(
                processors=console_processors,
                wrapper_class=structlog.make_filtering_bound_logger(
                    getattr(structlog.stdlib.logging, console_level_upper, structlog.stdlib.logging.INFO)
                ),
                context_class=dict,
                logger_factory=structlog.PrintLoggerFactory(file=sys.stdout),
                cache_logger_on_first_use=False,
            )

            self.logger = structlog.get_logger("validation_pipeline")

    def debug(self, message: str, **kwargs):
        """Log debug message with optional structured context."""
        if self.logger:
            self.logger.debug(message, **kwargs)

    def info(self, message: str, **kwargs):
        """Log info message with optional structured context."""
        if self.logger:
            self.logger.info(message, **kwargs)

    def warning(self, message: str, **kwargs):
        """Log warning message with optional structured context."""
        if self.logger:
            self.logger.warning(message, **kwargs)
        with self._issues_lock:
            self.validation_issues.append(('WARNING', message))

    def error(self, message: str, **kwargs):
        """Log error message with optional structured context."""
        if self.logger:
            self.logger.error(message, **kwargs)
        with self._issues_lock:
            self.validation_issues.append(('ERROR', message))

    def critical(self, message: str, **kwargs):
        """Log critical message with optional structured context."""
        if self.logger:
            self.logger.critical(message, **kwargs)
        with self._issues_lock:
            self.validation_issues.append(('CRITICAL', message))

    def add_validation_issue(self, level: str, category: str, message: str, details: dict = None):
        """Add a structured validation issue (thread-safe)."""
        issue = {
            'timestamp': datetime.now().isoformat(),
            'level': level,
            'category': category,
            'message': message,
            'details': details or {}
        }
        with self._issues_lock:
            self.validation_issues.append(issue)

        log_message = f"[{category}] {message}"
        log_kwargs = {'category': category}
        if details:
            log_kwargs.update(details)

        if self.logger:
            if level == 'ERROR':
                self.logger.error(log_message, **log_kwargs)
            elif level == 'WARNING':
                self.logger.warning(log_message, **log_kwargs)
            else:
                self.logger.info(log_message, **log_kwargs)

    def start_timer(self, name: str):
        """Start a named timer for performance measurement (thread-safe)."""
        with self._timers_lock:
            self._timers[name] = time.time()

    def stop_timer(self, name: str) -> float:
        """Stop a named timer and return elapsed time in seconds (thread-safe)."""
        end_time = time.time()
        with self._timers_lock:
            if name not in self._timers:
                raise KeyError(f"Timer '{name}' was never started")
            start_time = self._timers[name]
            elapsed = end_time - start_time
            self._timers[name] = elapsed
            return elapsed

    def get_timers(self) -> Dict[str, float]:
        """Get all recorded timers (thread-safe)."""
        with self._timers_lock:
            return self._timers.copy()

    def add_file_timing(self, input_file: str, validator_type: str, elapsed_time: float):
        """Add file timing information (thread-safe)."""
        timing = FileTimingSummary(
            input_file=input_file,
            validator_type=validator_type,
            elapsed_time=elapsed_time
        )
        with self._timings_lock:
            self.file_timings.append(timing)

    def clear_file_timings(self):
        """Clear all file timing information (thread-safe)."""
        with self._timings_lock:
            self.file_timings.clear()

    def display_file_timings_summary(self):
        """Display a formatted summary of file processing times."""
        if not self.file_timings:
            return

        with self._timings_lock:
            timings_copy = self.file_timings.copy()

        total_time = sum(t.elapsed_time for t in timings_copy)

        self.info("")
        self.info("=" * 80)
        self.info("FILE PROCESSING SUMMARY")
        self.info("=" * 80)

        max_filename_len = max(len(t.input_file) for t in timings_copy) if timings_copy else 20
        max_filename_len = min(max_filename_len, 50)

        for timing in timings_copy:
            filename = timing.input_file
            if len(filename) > max_filename_len:
                filename = "..." + filename[-(max_filename_len-3):]

            time_str = f"{timing.elapsed_time:.2f}s"
            type_label = timing.validator_type.upper()

            self.info(
                f"  {filename:<{max_filename_len}}  [{type_label:>7}]  {time_str:>8}",
                category=timing.validator_type
            )

        self.info("-" * 80)
        self.info(f"  {'TOTAL':<{max_filename_len}}             {total_time:>8.2f}s")
        self.info("=" * 80)
        self.info("")

    def clear_issues(self):
        """Clear all validation issues (thread-safe)."""
        with self._issues_lock:
            self.validation_issues.clear()

    def __enter__(self):
        """Context manager entry - clear issues at start."""
        self.clear_issues()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit - propagate exceptions."""
        return False


def setup_logging(
    console_level: str = "INFO",
    log_file: Optional[Path] = None,
) -> ValidationLogger:
    """Create and set up a new ValidationLogger instance."""
    logger = ValidationLogger()
    logger.setup(console_level, log_file)
    return logger


__all__ = ['ValidationLogger', 'setup_logging', 'get_incremented_path']
