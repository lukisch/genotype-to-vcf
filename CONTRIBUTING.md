# Contributing to Genotype-to-VCF Pro Converter

Thank you for your interest in contributing!

## How to Contribute

### Reporting Bugs

- Open a [GitHub Issue](../../issues) with a clear description
- Include your OS, Python version, and the error message
- **Never** include personal genetic data in bug reports

### Feature Requests

- Open a [GitHub Issue](../../issues) with the label `enhancement`
- Describe the use case and expected behavior

### Pull Requests

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/my-feature`)
3. Make your changes to `Make23toVCF3.py`
4. Test with sample data (not real genetic data)
5. Commit your changes (`git commit -m 'Add feature'`)
6. Push to the branch (`git push origin feature/my-feature`)
7. Open a Pull Request

### Code Style

- Python 3.8+ compatible
- Follow PEP 8 conventions
- Keep the single-file architecture (all logic in `Make23toVCF3.py`)
- Comment complex bioinformatics logic

### Privacy

- **Never** commit real genetic data (VCF, raw genotype files)
- Use synthetic test data for testing
- Do not include personal identifiers in code or commits

## Development Setup

```bash
pip install -r requirements.txt
python Make23toVCF3.py
```

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
