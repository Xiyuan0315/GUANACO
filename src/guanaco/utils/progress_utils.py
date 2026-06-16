"""Progress indicators for GUANACO startup"""

import sys
import threading
import time

class Spinner:
    """Animated spinner for long operations"""
    
    def __init__(self, message="Loading", *, show_result=True, stream=None):
        self.message = message
        self.running = False
        self.thread = None
        self.spinner_chars = ["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"]
        self.current = 0
        self.show_result = show_result
        self.stream = stream or sys.stdout
        
    def __enter__(self):
        self.start()
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.stop()
        
    def start(self):
        """Start the spinner"""
        self.running = True
        self.thread = threading.Thread(target=self._spin)
        self.thread.daemon = True
        self.thread.start()
        
    def _spin(self):
        """Spin animation"""
        while self.running:
            char = self.spinner_chars[self.current % len(self.spinner_chars)]
            self.stream.write(f'\r{char} {self.message}')
            self.stream.flush()
            self.current += 1
            time.sleep(0.1)
            
    def stop(self, success=True):
        """Stop the spinner"""
        self.running = False
        if self.thread:
            self.thread.join()
        
        # Clear the line and show result
        self.stream.write('\r' + ' ' * (len(self.message) + 4) + '\r')
        if self.show_result:
            if success:
                self.stream.write(f"✓ {self.message}\n")
            else:
                self.stream.write(f"✗ {self.message}\n")
        self.stream.flush()


def progress_bar(current, total, message="", width=40):
    """Simple progress bar for known quantities"""
    percent = current / total
    filled = int(width * percent)
    bar = "█" * filled + "░" * (width - filled)
    sys.stdout.write(f'\r{message} [{bar}] {percent*100:.0f}%')
    sys.stdout.flush()
    if current >= total:
        print()  # New line when complete
